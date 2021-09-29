# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 16:38:37 2021

@author: ashoff
"""

import os
import re
import glob2
import math
import errno
from shutil import copyfile


# Timestamps
from datetime import datetime as dt

# Data organization
import pandas as pd
import numpy as np

# Data Parsing
from itertools import islice

# Data Fitting
from lmfit import Parameters, minimize

# X-ray Science
import larch
import xraylib

# Atomic Environment Geneartion and Manipulation
from ase.visualize import view

# Plotting
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec


##############################################################################

                        # NON-CLASS GENERIC FUNCTIONS #
                        
##############################################################################


def calc_mu(numerator, denominator, log=True, flip = False):
    """
    Calculates the adsorption coeffieicnt mu from two arrays, numerator and
    denomiator. mu = N/D, ln(N/D), D/N, or ln(D/N) depending on args.

    Parameters
    ----------
    numerator : list or array
        array of size n to use as numerator
    denominator : list or array
        array of size n to use as denominator
    log : bool, optional
        if True, take ln() of N/D. The default is True.
    flip : bool, optional
        Flip numerator and denominator. The default is False.

    Returns
    -------
    mu : numpy array
        arry of size n

    """
        
    if not flip:
        if not log:
            mu = np.divide(numerator, denominator)
        elif log:
            mu = np.log(np.divide(numerator, denominator))
        else:
            print("ERROR: SET LOG TO BOOL")
        
    elif flip:
        if not log:
            mu = np.divide(numerator, denominator)
        elif log:
            mu = np.log(np.divide(numerator, denominator))
        else:
            print("ERROR: SET LOG TO BOOL")
        
    else:
        print("ERROR: SET FLIP TO BOOL")
            
    return mu

def calculate_spectrum_e0(larch_group, edge_energy, energy_range = 20):
        '''
        finds the edge positions of the samples and sets e0 to its value
        spectra = 'mu Sample'/'mu Reference'
        edge_energy = what energy should the edge be at
        energy_range = +/- range to look for inflection point

        returns a list of all edge energies found

        '''
        ind_min = find_nearest(larch_group.energy, edge_energy-20)[0]
        ind_max = find_nearest(larch_group.energy, edge_energy+20)[0]

        larch.xafs.find_e0(larch_group.energy[ind_min:ind_max], 
                           mu = larch_group.mu[ind_min:ind_max],
                           group = larch_group)

        e0_pos = larch_group.e0

        return e0_pos

def create_larch_spectrum(photon_energy, numerator, denominator, log=True, flip = False, name = None):
        '''
        docstring here, list, list, list, bool, bool
        '''
       
        spectrum = larch.Group()
        spectrum.energy = photon_energy
        spectrum.mu = calc_mu(numerator, denominator, log=True, flip = False)
        
        if name != None:
            spectrum.name = name
        
        return spectrum

def create_subdir(parent_dir, string):
    """
    Fuction to make a subdirectory in specified directory. Checks to see if 
    the directory already exists, and eitehr makes it or doesn't. Updates
    command line to the status of the subdirectory.
    
    Parameters
    ----------
    parent_dir : str
        directory to place subdirector
    string : str
        name of subdirectory

    Returns
    -------
    newdir : str
        full path strin of created director

    """

    newdir = os.path.join(parent_dir, string)

    try:
        os.mkdir(newdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            print("Directory Already Exists - Continuing Program")
        else:
            print ("Creation of the directory failed")
    else:
        print ("Successfully created subdirectory")
        
    return newdir
    
def CXAS_Sorted(files_directory, time_stamp = True, time_line = 0, time_format = '%m/%d/%Y %I:%M:%S %p', padded = True):
    """
    files: list of paths
    time [bool]: True if data contains headers that is a time stamp
    padded [bool]: true of data contained padded zeros, doesnt matter if time = True
    """
    
    files = glob2.glob(files_directory + '/*')
    
    path_series = pd.Series(files)
    
    filename_list = []
    padded_list = []
    time_list = []
    
    for line in files:
        
        filename = os.path.basename(line)[:-4]
        filename_list.append(filename)
        
        if not padded:
            scan_no = get_trailing_number(filename)
            char = len(str(scan_no))
            filename = filename[:-char] + str(scan_no).zfill(4)
            
        padded_list.append(filename)
        
        if time_stamp:
            # Open the file in read only mode
            with open(line, 'r') as f:
                count = 0
                # Read all lines in the file one by one
                for line in f:
                    # For each line, check if line contains the string
                    if count == time_line:
                        #print(line)
                        date = dt.strptime(line[:-1], time_format)
                        #print(date)
                        time_list.append(date)
                        break
                    else:
                        count = count + 1



    filename_series = pd.Series(filename_list)
    padded_series = pd.Series(padded_list)
    
    
    if time_stamp:
        time_series = pd.Series(time_list)
        
        temp_dict = {'Time': time_series, 'File Name': filename_series, 'Padded Name': padded_series, 'Path': path_series}
        
        TOS = pd.DataFrame(temp_dict)
        TOS.sort_values(by = 'Time', ignore_index = True, inplace = True)
        TOS.set_index('Time', inplace = True)
    
    else:
        temp_dict = {'File Name': filename_series, 'Padded Name': padded_series, 'Path': path_series}
        
        TOS = pd.DataFrame(temp_dict)
        TOS.sort_values(by = 'Padded Name', ignore_index = True, inplace = True)
        TOS.index.rename('Scan', inplace = True)
        
    return TOS

def find_nearest(array, value): 
    """
    Finds the index of the array that is closest to the requested value
    
    Returns: [index, array[index]]
    """
    array = np.asarray(array)
    
    idx = (np.abs(array - value)).argmin()
    
    return [idx, array[idx]]

def get_trailing_number(s):
    """
    Checks to see if there is a number at the end of string
    
    Parameters
    ----------
    s : string
            String to look for numbers at the end of

    Returns
    -------
    int
            Rumber at end of string, if none, return None

    """
    m = re.search(r'\d+$', s)
    return int(m.group()) if m else None
    
def interp_df(df, new_index):
    '''
    NEEDS Updatinng - index must be energy scale

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    new_index : TYPE
        DESCRIPTION.

    Returns
    -------
    df_out : TYPE
        DESCRIPTION.

    '''
    """Return a new DataFrame with all columns values interpolated
    to the new_index values."""
    df_out = pd.DataFrame(index=new_index)
    df_out.index.name = df.index.name

    for colname, col in df.iteritems():
        df_out[colname] = np.interp(new_index, df.index, col)

    return df_out

def interpolate_spectrum(larch_group, new_energy_axis):
    '''
    details here
    '''
    temp_dict = {'mu':larch_group.mu}
    index = larch_group.energy
    
    temp_df = pd.DataFrame(data = temp_dict, index = index)
    temp_df.index.rename('Energy', inplace = True)
    
    interpEnergy_df = interp_df(temp_df, new_energy_axis)
    
    larch_group.energy = new_energy_axis
    larch_group.mu = interpEnergy_df['mu'].values
    
    temp_dict.clear()
    
    return

def mergeindex(df1, df2, method = 'time'):
    """
    interpolates df2 onto df1 index
    
    df1: dataframe wiht index used for reindexing. requires datetime index
    
    df2: dataframe to be reindexed. requires datetime index
    
    Note: if df1 index exceeds limits of df2, NaN will be returned in values that can not be interpoalted
    
    Return: df2 wiht df1 index 
    """
    
    
    df2 = df2.reindex(df2.index.union(df1.index)).interpolate(method='time').reindex(df1.index)
    
    return df2
    
def ReadLVData(filename):
    """
    Reads a text file from the Co-ACCESS ambient pressure (AmP) labview system
    
    filename: string to file containing AmP labView data
    
    Returns: Pandas dataframe with a datetime index
    """
    # Extract Data into pandas dataframe
    data = pd.read_csv(filename, sep = '\t', encoding= 'unicode_escape')
    data.drop(data.columns[data.columns.str.contains('unnamed',case = False)],
              axis = 1, inplace = True)
    
    # Convert "Date and Time" to timedelta + date and time file was created
    data.rename(columns={'Date and Time':'Time'}, inplace=True)
    data['Time'] = pd.to_datetime(data['Time'])
    
    # Reindex the data to the time stamp     
    data = data.set_index('Time')
    
    return data

def ReadMSData(filename):
    """    
    Reads a CSV file supplied from the  Hiden Massoft figure file export.
    
    filename: string to file containing mass spectrometer data
    
    Returns: Pandas dataframe with a datetime index and signals matched with signal value (m/z) 
    """
    
    #Characteristics of Hiden Data Structure
    header_line = 1
    date_line = 2
    scans_line = 4

    # Get Metadata Header info and store it in list
    with open(filename, 'r') as myfile:
        for i, line in enumerate(myfile):
            if i == header_line:
                header = [int(word) for word in line.rstrip().split(',') if word.isdigit()][0]
            elif i == date_line:
                date_str = [str(word) for word in line.rstrip().split(',')]
                date = dt.strptime(date_str[1] + ' ' + date_str[3],'%m/%d/%Y %I:%M:%S %p')
            elif i == scans_line:
                scans = [int(word) for word in line.rstrip().split(',') if word.isdigit()][0]
            elif i > 50:
                break

    myfile.close()
    
    # Get Details about numebr of scans and their ranges
    Scan_List = [['Scan Number', 'Input', 'Start', 'Stop', 'Increment']]
    with open(filename, 'r') as myfile:
        for line in islice(myfile, scans_line + 2, scans_line + 2 + scans):
            temp_scan = line.strip().split(',')
            Scan_List.append([temp_scan[0], temp_scan[2], temp_scan[4], temp_scan[5], temp_scan[6]])
        
    myfile.close()
    
    '''# Feedback the key header characteristics
    print('Header Lines: {0}'.format(header))
    print('Scan Start Time: {0}'.format(date))
    print('Number of scans: {0}'.format(scans))
    for line in Scan_List:
        print(line)'''    
    
    # Extract Data into pandas dataframe
    data = pd.read_csv(filename, sep = ',', header = header+1)
    data.drop(data.columns[data.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
    
    # Convert timestamp to timedelta + date and time file was created
    data['Time'] = pd.to_timedelta(data['Time']) + date
    data.drop('ms', axis = 1, inplace = True)
    
    '''# Feedback to see column names
    print('Column Headers = {0}'.format(list(data.columns)))'''
    
    
    # Conver Scan_List parameters into list of headers to apply to data    
    if all(x in list(data.columns) for x in ['Cycle', 'mass amu']): # Implies there is 'Cycle' Data and is processed accoridngly
        # Defines Column Headers
        min_amu = int(float(Scan_List[1][2]))
        max_amu = int(float(Scan_List[1][3]))
        col_names = list(range(min_amu, max_amu+1))
        
        '''# Feedback to see that it gets the AMU ranges correct:
        print('Start amu value = {0}'.format(min_amu))
        print('End amu value = {0}'.format(max_amu))
        print('New Column Lists = {0}'.format(col_names))'''
        
        # Restructures dataframe into select m/z values as f(time)
        data.drop(['Cycle'], axis = 1, inplace = True)
        data.set_index(['mass amu','Time'], inplace = True)
        for i, line in enumerate(col_names):
            if i == 0:
                df = data.xs(line)['SEM Torr']
            else:
                temp_df = data.xs(line)['SEM Torr']
                df = pd.concat([df, temp_df], axis = 1)
        
        # Relables Columns
        df.set_axis(col_names, axis = 1, inplace = True)
        data = df
          
    else:
        if "Cycle" in list(data):
            data.drop('Cycle', axis = 1, inplace = True)
        
        col_names = list(data)
        
        new_col_names = []
        for line in col_names:
            string = line.strip()
            if re.match(r'^-?\d+(?:\.\d+)$', string):
                new_col_names.append(float(string))
            else:
                new_col_names.append(line)
        data.columns = new_col_names
        data.set_index('Time', inplace = True)
    
    return data



##############################################################################

                        # NON-CLASS LCF FUNCTIONS #
                        
##############################################################################

def sum_standards(pars, basis_spectra):
    '''
    pars is paramters group with number of elemnts equal to the number of basis
    basis is a df with rows = energy, columns = basis#, column names = basis#
    data is larch group of same lenght of basis
    '''    
    #consider making a basis set a df not a larh group
    model = np.zeros(len(basis_spectra.index))

   
    for amp, basis in zip(pars.keys(), list(basis_spectra.columns)):
        model = model + pars[amp].value*basis_spectra[basis].values
    
    return model

# Define function to calcualte residuals
def resid(pars, data, basis):
    '''
    pars is paramters group
    data is a list
    basis is a df with rows = energy, columns = basis#, column names = basis#
    '''
    
    #consider simplifying for lsits or dfs not larch groups
    residual = (data - sum_standards(pars, basis))#/data.eps
    
    return residual


##############################################################################

            # EXPERIMENTAL CLASS - IN-SITU XAS PROCESSING #
                        
##############################################################################

class Experiment:
    # Constructor
    def __init__(self, name):
        '''
        docstring here
        '''
        self.name = name
        # Place to store all process related data streams
        self.process_params = {} 
        #P lace to store all beamline data and XAS workup
        self.spectra = {} 
        # Place to store all  analysis results from things like LCF/PCA
        self.analysis = {'LCF':{}, 'PCA':{}} 
         # Catch-all place for summary of results, needs work
        self.summary = {}
        
        return
    
    # Functions
    def import_massspec(self, file_name):
        '''
        docstring here
        '''
        self.process_params['MS Data'] = ReadMSData(file_name)
        self.summary['MS Filename'] = file_name
        
        return
        
    def import_labview(self, file_name):
        '''
        docstring here
        '''
        self.process_params['LV Data'] = ReadLVData(file_name)
        self.summary['LV Filename'] = file_name
        
        return
        
    def import_spectra_data(self, xas_data_directory, xas_data_structure):
        '''
        docstring here, list, dictionary
        '''
        self.summary['XAS Data Structure'] = xas_data_structure
        
        time_stamp = self.summary['XAS Data Structure']['time stamp']
        time_line = self.summary['XAS Data Structure']['time on line']
        time_format = self.summary['XAS Data Structure']['time format']
        padded = self.summary['XAS Data Structure']['padded scan numbers']
        
        self.summary['XAS Spectra Files'] = CXAS_Sorted(xas_data_directory, time_stamp = time_stamp, time_line = time_line, time_format = time_format, padded = padded)  
        
        sep = self.summary['XAS Data Structure']['separator']
        names = self.summary['XAS Data Structure']['column names']
        skiprows = self.summary['XAS Data Structure']['skiprows']
        eneryg_index = self.summary['XAS Data Structure']['energy column']
        
        
        for index, row in self.summary['XAS Spectra Files'].iterrows():
            
            filename = row[0]
            file_path = row[2]

            self.spectra[filename] = {}
        
            self.spectra[filename]['Time'] = index
            
            self.spectra[filename]['BL Data'] = pd.read_csv(file_path, names = names, sep = sep, skiprows = skiprows, index_col = False)
        
        return
     
    def correlate_process_params(self):
        '''
        update me here, looks through every dictionary in process params and tries to interpolate the parameter values
        onto a spectra time stamp
        '''          
        
        self.summary['XAS Spectra Process Params'] = self.summary['XAS Spectra Files']['File Name']
        
        for key in self.process_params.keys():                
            temp_df = mergeindex(self.summary['XAS Spectra Process Params'], self.process_params[key])
            
            self.summary['XAS Spectra Process Params'] = pd.concat([self.summary['XAS Spectra Process Params'], temp_df], axis=1)
            
        for key in self.spectra.keys():
            self.spectra[key]['Process Values'] =  self.summary['XAS Spectra Process Params'].loc[[self.spectra[key]['Time']]]
        
        return
    
    def calculate_spectra(self):
        '''
        docstring TBD
        '''
        energy_column = self.summary['XAS Data Structure']['energy column']
        sample_numerator = self.summary['XAS Data Structure']['sample numerator']
        sample_denominator = self.summary['XAS Data Structure']['sample denominator']
        sample_ln = self.summary['XAS Data Structure']['sample ln']
        sample_invert = self.summary['XAS Data Structure']['sample invert'] 
        
        reference_numerator = self.summary['XAS Data Structure']['reference numerator']
        reference_denominator = self.summary['XAS Data Structure']['reference denominator']
        reference_ln = self.summary['XAS Data Structure']['reference ln']
        reference_invert = self.summary['XAS Data Structure']['reference invert'] 
        
        for key in self.spectra.keys():
            photon_energy = self.spectra[key]['BL Data'][energy_column].values
            samp_numerator = self.spectra[key]['BL Data'][sample_numerator].values
            samp_denominator = self.spectra[key]['BL Data'][sample_denominator].values
            samp_log=sample_ln
            samp_flip = sample_invert
            ref_numerator = self.spectra[key]['BL Data'][reference_numerator].values
            ref_denominator = self.spectra[key]['BL Data'][reference_denominator].values
            ref_log = reference_ln
            ref_flip = reference_invert

            self.spectra[key]['mu Sample'] = create_larch_spectrum(photon_energy, samp_numerator, samp_denominator, log=samp_log, flip = samp_flip)
            self.spectra[key]['mu Reference'] = create_larch_spectrum(photon_energy, ref_numerator, ref_denominator, log=ref_log, flip = ref_flip)
        return    

    def interpolate_spectra(self, offset = 5):
        '''
        NEEDS UPDATING 

        Parameters
        ----------
        spectra_dict : TYPE
            DESCRIPTION.
        offset : TYPE, optional
            DESCRIPTION. The default is 5.

        Returns
        -------
        None.

        '''
        # Empty lists to store energy values
        low_limit = []
        high_limit = []

        # Find the highest low and lowest high for all energy scales in the dictionary
        for key in self.spectra.keys():
            low_limit.append(self.spectra[key]['mu Sample'].energy.min())
            high_limit.append(self.spectra[key]['mu Sample'].energy.max())

        print('Starting Energy Statistics:')
        print(f'\t Minimim Starting Energy: {min(low_limit):.2f}')
        print(f'\t Maximum Starting Energy: {max(low_limit):.2f}')
        print(f'\t Starting Energy mean: {np.asarray(low_limit).mean():.2f} +/- {np.asarray(low_limit).std():.2f}')

        print('\nEnding Energy Statistics:')
        print(f'\t Minimim Ending Energy: {min(high_limit):.2f}')
        print(f'\t Maximum Ending Energy: {max(high_limit):.2f}')
        print(f'\t Ending Energy mean: {np.asarray(high_limit).mean():.2f} +/- {np.asarray(high_limit).std():.2f}')

        # Find index of first spectra that are bounded by the offset
        offset = offset
        temp_key = next(iter(self.spectra)) 
        min_ind = find_nearest(self.spectra[temp_key]['mu Sample'].energy, max(low_limit)+offset)[0]
        max_ind = find_nearest(self.spectra[temp_key]['mu Sample'].energy, min(high_limit)-offset)[0]

        print(f'\nEnergy Range Interpolated Between: {self.spectra[temp_key]["mu Sample"].energy[min_ind]:.2f}, {self.spectra[temp_key]["mu Sample"].energy[max_ind]:.2f}')
        print(f'Number of Data Points: {max_ind-min_ind+1}')

        # Determine the new index to interpoalte onto from the index limits above
        CommonEnergy_index = self.spectra[temp_key]['mu Sample'].energy[min_ind:max_ind]

        #Interpolate all data sets onto new axis, and build an interpoalted dictionary
        for key in self.spectra.keys():
            interpolate_spectrum(self.spectra[key]['mu Sample'], CommonEnergy_index)
            interpolate_spectrum(self.spectra[key]['mu Reference'], CommonEnergy_index)
        
        return

    
    def calibrate_reference_spectra(self, edge_energy, energy_range=20):
        '''
        Calibrates the reference channel spectra based upon edge_energy provided
        sets reference energy e0 to edge_energy
        shifts sample + reference energy scales based upon average edge energy position
        '''
        e0_list = []
        
        for key in self.spectra.keys():
            e0_list.append(calculate_spectrum_e0(self.spectra[key]['mu Reference'], edge_energy, energy_range = energy_range))        

        #Shift energy scale based upon E0 values
        del_E = edge_energy - np.asarray(e0_list).mean()

        for key in self.spectra.keys():
            self.spectra[key]['mu Reference'].e0 = edge_energy
            self.spectra[key]['mu Reference'].energy = self.spectra[key]['mu Reference'].energy + del_E
            self.spectra[key]['mu Sample'].energy = self.spectra[key]['mu Sample'].energy + del_E
        
        print('Reference Calibraiton Statistics:')
        print(f'Reference E0 min: {min(e0_list):.2f} eV')
        print(f'Reference E0 max: {max(e0_list):.2f} eV')
        print(f'Reference E0 mean: {np.asarray(e0_list).mean():.2f} +/- {np.asarray(e0_list).std():.2f} eV')
        print(f'Reference E0 calibrated to: {edge_energy:.2f} eV')
        print(f'Spectra shifted by {del_E:.2f} eV\n\n')

        return
        
    def find_sample_e0(self, edge_energy, energy_range = 20, use_mean = True):
        '''
        finds the edge positions of the samples around edge_energy value and energy_range
        set sample e0 value to found value or mean value
        use_mean = True --> use mean value of e0 for every spectra
        use_mean = False --> use calculated value of e0 for each spectra

        '''

        e0_list = []
        
        for key in self.spectra.keys():
            e0_list.append(calculate_spectrum_e0(self.spectra[key]['mu Sample'], edge_energy, energy_range = energy_range))

        if use_mean:
            for key in self.spectra.keys():
                self.spectra[key]['mu Sample'].e0 = np.asarray(e0_list).mean()
        
        print('Sample Calibraiton Statistics:')
        print(f'Sample E0 min: {min(e0_list):.2f} eV')
        print(f'Sample E0 max: {max(e0_list):.2f} eV')
        print(f'Sample E0 mean: {np.asarray(e0_list).mean():.2f} +/- {np.asarray(e0_list).std():.2f} eV')
        if use_mean:
            print(f'Sample E0 set to: {np.asarray(e0_list).mean():.2f} eV\n\n')
        else:
            print('Sample E0 calculated for each spectra\n\n')
        
        return
    
    def set_normalization_parameters(self, spectra_name, pre1 = -100, pre2 = -50, norm1 = 75, norm2 = 300, nnorm = 2, make_flat = True):
        
        for key in self.spectra.keys():
            self.spectra[key][spectra_name].pre1 = pre1
            self.spectra[key][spectra_name].pre2 = pre2
            self.spectra[key][spectra_name].norm1 = norm1
            self.spectra[key][spectra_name].norm2 = norm2
            self.spectra[key][spectra_name].nnorm = nnorm
            self.spectra[key][spectra_name].make_flat = make_flat
        
        return
        
    def normalize_spectra(self, spectra_name):        
        '''
        FILL ME IN
        spectra

        Returns
        -------
        None.

        '''
                
        for key in self.spectra.keys():
            energy = self.spectra[key][spectra_name].energy
            mu = self.spectra[key][spectra_name].mu
            group = self.spectra[key][spectra_name]
            e0 = self.spectra[key][spectra_name].e0
            pre1 = self.spectra[key][spectra_name].pre1
            pre2 = self.spectra[key][spectra_name].pre2
            norm1 = self.spectra[key][spectra_name].norm1
            norm2 = self.spectra[key][spectra_name].norm2
            nnorm = self.spectra[key][spectra_name].nnorm
            make_flat = self.spectra[key][spectra_name].make_flat
            
            larch.xafs.pre_edge(energy, mu = mu, group = group, e0 = e0, 
                                pre1 = pre1, pre2 = pre2, 
                                norm1 = norm1, norm2 = norm2, 
                                nnorm = nnorm, make_flat = make_flat)

        return
    
    
       
    ### LCF related functions
    def load_lcf_basis(self, basis_list):
        '''
        update later basis_must be a list of larch groups wiht nornalized spectra
        '''
        
        # Loads in the basis spectra
        self.analysis['LCF']['basis spectra'] = {}   
        
        for i in range(len(basis_list)):
            self.analysis['LCF']['basis spectra']['basis'+  str(i +1)] = basis_list[i]
        
        # Creates a paramters space for each spectra in self.spectra
        total_spectra = len(self.spectra.keys()) 
        energy_range = total_spectra*[None]
        parameters = total_spectra*[None]
        fitting_parameters = total_spectra*[None]
        names = list(self.spectra.keys())
        
        self.analysis['LCF']['fitting parameters'] = pd.DataFrame({'names': names,'energy range':energy_range, 'lcf parameters':parameters, 'lcf report':fitting_parameters})
        self.analysis['LCF']['fitting parameters'].set_index('names', inplace = True)
        
        for ind in self.analysis['LCF']['fitting parameters'].index:
            
            self.analysis['LCF']['fitting parameters']['lcf parameters'][ind] = Parameters()
            
            for x in range(len(self.analysis['LCF']['basis spectra'].keys())):
                 
                 self.analysis['LCF']['fitting parameters']['lcf parameters'][ind].add(f'amp{x+1}', value = 1/len(basis_list), min = 0, max = 1, vary = True)
        
        return    
    
    def fit_LCF(self, energy_range):
        '''
        energy_range = [emin, emax]
        '''
        for ind in self.analysis['LCF']['fitting parameters'].index:
            self.analysis['LCF']['fitting parameters']['energy range'][ind] = energy_range

        for key in self.spectra.keys():

            i1 = find_nearest(self.spectra[key]['mu Sample'].energy, energy_range[0])[0]
            i2 = find_nearest(self.spectra[key]['mu Sample'].energy, energy_range[1])[0]

            data_to_fit = self.spectra[key]['mu Sample'].flat[i1:i2+1]
            
            basis_dict = {}


            for key2 in self.analysis['LCF']['basis spectra'].keys():
                basis_dict[key2] = self.analysis['LCF']['basis spectra'][key2].flat[i1:i2+1]

            basis_data = pd.DataFrame(basis_dict)
            

            result = minimize(resid, self.analysis['LCF']['fitting parameters']['lcf parameters'][key],
                              args=(data_to_fit, basis_data,))

            self.analysis['LCF']['fitting parameters']['lcf report'][key] = result
            

            # Add fit summary to each spectra file 'LCF Results'
            energy = self.spectra[key]['mu Sample'].energy
            spectra = self.spectra[key]['mu Sample'].flat
            
            LCF_result_dict = {'Energy': energy, 'Spectra': spectra}
            
            basis_data_dict = {}
            
            for amp, key2 in zip(self.analysis['LCF']['fitting parameters']['lcf report'][key].params, self.analysis['LCF']['basis spectra'].keys()):
                weighting_fraction = self.analysis['LCF']['fitting parameters']['lcf report'][key].params[amp].value
                basis_name = self.analysis['LCF']['basis spectra'][key2].name
                basis_data_dict[key2 + '_'+ basis_name] = self.analysis['LCF']['basis spectra'][key2].flat
                LCF_result_dict[key2 + '_'+ basis_name] = weighting_fraction*self.analysis['LCF']['basis spectra'][key2].flat

            basis_data = pd.DataFrame(basis_data_dict) 
        
            fit = sum_standards(self.analysis['LCF']['fitting parameters']['lcf report'][key].params, basis_data)
            LCF_result_dict['Fit'] = fit
        
            residuals = resid(self.analysis['LCF']['fitting parameters']['lcf report'][key].params,
                  self.spectra[key]['mu Sample'].flat, basis_data)
        
            LCF_result_dict['Residuals'] = residuals

            lcf_result_df = pd.DataFrame(LCF_result_dict)

            lcf_result_df.set_index('Energy', inplace = True)
            
            self.spectra[key]['LCF Results'] = lcf_result_df
        
        return
   
    def lcf_report(self):
    
        no_basis = len(self.analysis['LCF']['basis spectra'].keys())

        fit_param = {'Name': [],
                    'Chi2': [],
                    'RedChi2': [],
                    'Variables': []}

        for x in range(no_basis):
            fit_param[f'Amp{x+1}'] = []
            fit_param[f'Amp{x+1}-stdev'] = []
            
        fit_param['Sum Amp'] = []
            

        for ind in self.analysis['LCF']['fitting parameters']['lcf report'].index:
            fit_param['Name'].append(ind)
            fit_param['Chi2'].append(self.analysis['LCF']['fitting parameters']['lcf report'][ind].chisqr)
            fit_param['RedChi2'].append(self.analysis['LCF']['fitting parameters']['lcf report'][ind].redchi)
            fit_param['Variables'].append(self.analysis['LCF']['fitting parameters']['lcf report'][ind].nvarys)

                    
            sum_amp = 0
            
            for x in range(no_basis):
                fit_param[f'Amp{x+1}'].append(self.analysis['LCF']['fitting parameters']['lcf report'][ind].params[f'amp{x+1}'].value)
                fit_param[f'Amp{x+1}-stdev'].append(self.analysis['LCF']['fitting parameters']['lcf report'][ind].params[f'amp{x+1}'].stderr)
                sum_amp = sum_amp + self.analysis['LCF']['fitting parameters']['lcf report'][ind].params[f'amp{x+1}'].value
            
            fit_param['Sum Amp'].append(sum_amp)
                
        LCF_df = pd.DataFrame(fit_param)

        self.analysis['LCF']['Fit Summary'] = LCF_df
    
    
        return

    ### Data Exporting Functions
    def save_normalize_spectra(self, output_directory):
        # Define filename of results:
        fname_normspectra = self.name + '_NormalizedSpectra.csv'
        
        # Define where the data will be saved
        output_path = os.path.join(output_directory, fname_normspectra)
        
        for i, key in zip(range(len(self.spectra.keys())), list(self.spectra.keys())):
            if i == 0:
                normalized_df = pd.DataFrame({'Energy': self.spectra[key]['mu Sample'].energy, f'{key}': self.spectra[key]['mu Sample'].flat})
                normalized_df.set_index('Energy', inplace = True)
            else:
                temp_df = pd.DataFrame({'Energy': self.spectra[key]['mu Sample'].energy, f'{key}': self.spectra[key]['mu Sample'].flat})
                temp_df.set_index('Energy', inplace = True)
                
                normalized_df = pd.concat([normalized_df, temp_df], axis = 1)
        
        normalized_df.to_csv(output_path, sep=',', na_rep='', header=True, index=True)
        
        return
        
        
    
    def save_processparams(self, output_directory):
        # Define filename of results:
        fname_processparams = self.name + '_ProcessPrams.csv'

        # Define where the data will be saved
        output_path = os.path.join(output_directory, fname_processparams)
    
        # Save the data
        self.summary['XAS Spectra Process Params'].to_csv(output_path, sep=',', na_rep='', header=True, index=True)
    
        print('Process Parameter Data Saved')
    
        return
    
    def save_lcf_results(self, output_directory, save_spectrum = True):
        # Define filename of results:
        fname_LCFSummary = self.name + '_LCFsummary.csv'
    
        # Define where the data will be saved
        output_path = os.path.join(output_directory, fname_LCFSummary)
    
        # Save the LCF summary data
        self.analysis['LCF']['Fit Summary'].to_csv(output_path, sep=',', na_rep='', header=True, index=True)
    
        if save_spectrum:    
            # Create Subdirectory to store individual fit results
            sub_directory = create_subdir(output_directory, 'Spectra LCF Results')
        
            
            for key in self.spectra.keys():
                # Define filename of results:
                fname_spectrumLCF = key + '_LCFsummary.csv'
                
                # Define where the data will be saved
                output_path2 = os.path.join(sub_directory, fname_spectrumLCF)
                
                # Save the LCF summary data
                self.spectra[key]['LCF Results'].to_csv(output_path2, sep=',', na_rep='', header=True, index=True)
                
        print('LCF Data Saved')
    
        return
        
    
    
    ### Plotting Data - needs improvement
    def plot_xas_spectra(self, emin, emax, calib_line = False):
        # Define Figure [2 panel side by side]
        fig1 = plt.figure(constrained_layout=True, figsize = (12,5))
        spec1 = gridspec.GridSpec(ncols = 2, nrows = 1, figure = fig1)

        f1_ax1 = fig1.add_subplot(spec1[0])
        f1_ax2 = fig1.add_subplot(spec1[1])

        #Plot Reference and Sample Spectra
        for key in self.spectra.keys():
            f1_ax1.plot(self.spectra[key]['mu Reference'].energy, self.spectra[key]['mu Reference'].mu)
            f1_ax2.plot(self.spectra[key]['mu Sample'].energy, self.spectra[key]['mu Sample'].mu)
        
        if calib_line:
            f1_ax1.plot([self.spectra[key]['mu Reference'].e0, 
                         self.spectra[key]['mu Reference'].e0],
                        [min(self.spectra[key]['mu Reference'].mu), max(self.spectra[key]['mu Reference'].mu)],
                        color = 'k')
        emin_ref = emin
        emax_ref = emax

        emin_samp = emin
        emax_samp = emax

        f1_ax1.set_xlim([emin_ref, emax_ref])
        f1_ax1.set_title('Reference')
        f1_ax1.set_xlabel('Photon Energy (eV)')
        f1_ax1.set_ylabel('mu(E)x')

        f1_ax2.set_xlim([emin_samp, emax_samp])
        f1_ax2.set_title('Sample')
        f1_ax2.set_xlabel('Photon Energy (eV)')
        f1_ax2.set_ylabel('mu(E)x')
        
    def plot_norm_spectra(self):
    # Define Figure [2 panel side by side]
        fig1 = plt.figure(constrained_layout=True, figsize = (12,10))
        spec1 = gridspec.GridSpec(ncols = 2, nrows = 2, figure = fig1)

        f1_ax1 = fig1.add_subplot(spec1[0,0])
        f1_ax2 = fig1.add_subplot(spec1[0,1])
        f1_ax3 = fig1.add_subplot(spec1[1,0])
        f1_ax4 = fig1.add_subplot(spec1[1,1])

        #Plot Reference and Sample Spectra
        for key in self.spectra.keys():
            f1_ax1.plot(self.spectra[key]['mu Reference'].energy, self.spectra[key]['mu Reference'].flat)
            f1_ax2.plot(self.spectra[key]['mu Sample'].energy, self.spectra[key]['mu Sample'].flat)
            f1_ax3.plot(self.spectra[key]['mu Reference'].energy, self.spectra[key]['mu Reference'].flat)
            f1_ax4.plot(self.spectra[key]['mu Sample'].energy, self.spectra[key]['mu Sample'].flat)


        emin_ref1 = self.spectra[key]['mu Reference'].e0 - 200
        emax_ref1 = self.spectra[key]['mu Reference'].e0 + 1000
        emin_ref2 = self.spectra[key]['mu Reference'].e0 - 50
        emax_ref2 = self.spectra[key]['mu Reference'].e0 + 120

        emin_samp1 = self.spectra[key]['mu Sample'].e0 - 200
        emax_samp1 = self.spectra[key]['mu Sample'].e0 + 1000
        emin_samp2 = self.spectra[key]['mu Sample'].e0 - 50
        emax_samp2 = self.spectra[key]['mu Sample'].e0 + 120

        f1_ax1.set_xlim([emin_ref1, emax_ref1])
        f1_ax1.set_title('Reference')
        f1_ax1.set_xlabel('Photon Energy (eV)')
        f1_ax1.set_ylabel('Norm. mu(E)x')

        f1_ax2.set_xlim([emin_samp1, emax_samp1])
        f1_ax2.set_title('Sample')
        f1_ax2.set_xlabel('Photon Energy (eV)')
        f1_ax2.set_ylabel('Norm. mu(E)x')


        f1_ax3.set_xlim([emin_ref2, emax_ref2])
        f1_ax3.set_title('Reference')
        f1_ax3.set_xlabel('Photon Energy (eV)')
        f1_ax3.set_ylabel('Norm. mu(E)x')


        f1_ax4.set_xlim([emin_samp2, emax_samp2])
        f1_ax4.set_title('Sample')
        f1_ax4.set_xlabel('Photon Energy (eV)')
        f1_ax4.set_ylabel('Norm. mu(E)x')
        
    def plot_LCF_results(self, process_parameter = None):
        # Define Figure [2 panel side by side]
        fig1 = plt.figure(constrained_layout=True, figsize = (10,10))
        spec1 = gridspec.GridSpec(ncols = 6, nrows = 10, figure = fig1)

        f1_ax1 = fig1.add_subplot(spec1[0:7,0:])
        f1_ax2 = fig1.add_subplot(spec1[7:,0:], sharex=f1_ax1)

        for col in self.analysis['LCF']['Fit Summary'].columns:
            if not 'stdev' in col:
                if 'Amp' in col:
                    f1_ax1.plot(self.analysis['LCF']['Fit Summary'][col], label = col)

        f1_ax2.plot(self.analysis['LCF']['Fit Summary']['RedChi2'], color = 'k')


        if process_parameter != None:
            f1_ax3 = f1_ax1.twinx()
            f1_ax3.plot(self.summary['XAS Spectra Process Params'][process_parameter].values, color = 'k')

        fig1.subplots_adjust(wspace=0, hspace=0)
        
        f1_ax1.legend(loc='center right')
        f1_ax1.set_ylabel('Fraction of Basis')
        f1_ax2.set_xlabel('Spectra')
        f1_ax2.set_ylabel('Reduced Chi^2')
        








##############################################################################

            # EXPERIMENTAL FUNCTIONS - STILL IN DEVELOPMENT #
                        
##############################################################################

def atoms2cluster(atoms, absorbing_atom_index, distance_cutoff, is_DFT = False, DFT_rep = [2,2,1]):
    '''
    This atoms2cluster function creates a bigger cluster from a cif input file
    used in conjunction wiht cluster2feffinp

    Parameters
    ----------
    atoms : str
        atoms object in cif file format
    absorbing_atom_index : str
        index of atom in the atoms input that is the "core"
    distance_cutoff : float
        longest distance that an atom show be from the center of the cluster

    Returns
    -------
    atoms_cluster : ASE Atoms
        Atoms type of all atoms in the cluster defined by distance_cutoff.
    absorbing_atom_index : int
        index of the absorbing [core] atom in the atoms_cluster

    ''' 
    if is_DFT == True:
        atoms2 = atoms.repeat(DFT_rep)
        cen_atom = atoms2[absorbing_atom_index].position
        print('haha it passed teh DFT check')
    else:    
        # repeat the unitl cell to the nearest even repat such that the the unit cell lenghts are at least 2x distance_cutoff
        rep = max([math.ceil(distance_cutoff/atoms.cell[0][0])*2, 
               math.ceil(distance_cutoff/atoms.cell[1][1])*2, 
               math.ceil(distance_cutoff/atoms.cell[2][2])*2])
        
        print(f"Crystal repated {rep} time")
    
        # expand the unit cell
        atoms2 = atoms.repeat(rep)
        #view(atoms2)
    
        # define the center of the expanded unit cell
        cen_rep = [rep/2, rep/2, rep/2]
        #print(cen_rep)
    
        # define the center atom position of the unit cell
        cen_atom = np.sum(np.multiply(atoms.cell, cen_rep), axis = 0) + atoms[absorbing_atom_index].position
        print(f"Absorbing Atom Position: {cen_atom}")
    
    # Find new absorbing atom index
    absorbing_atom_index= [a.index for a in atoms2 if np.array_equal(a.position, cen_atom)]
    #print(f"Index of Absorbing Atom: {absorbing_atom_index}")

    # Remove atoms from model that are at a distance greater than distance_cutoff from the indexed atom
    atoms_cluster = atoms2[[a.index for a in atoms2 if atoms2.get_distance(absorbing_atom_index, a.index, mic=True) < distance_cutoff]]
    #view(atoms_cluster)
    
    
    # Update new absorbing atom index
    absorbing_atom_index= [a.index for a in atoms_cluster if np.array_equal(a.position, cen_atom)]
    print(f"Index of Absorbing Atom: {absorbing_atom_index}")
    
    
    return atoms_cluster, absorbing_atom_index[0]


def Calibrate_Spectra(spectra_dict, Ref_E0_calib, Samp_E0_calib):
    '''
    Calibrates the reference spectra defined in a spectra_dict, generated by ReadCXAS
    Updates the dictionaty to include edge position information

    Parameters
    ----------
    spectra_dict : Dictionary 
        generated from ReadCXAS.
    Ref_E0_calib : int
        energy (eV) that defines the reference channel edge.
    Samp_E0_calib : int
        approximate (+/- 5eV) edge of the sample, best approximation is the tabulated edge energy for the element of interest.

    Returns
    -------
    None.

    '''
    # Placeholder to determein statistics of reference edge energies
    Ref_E0_pos = []

    # Find Edge Energy of each Reference Spectra
    for key in spectra_dict.keys():
        ind_min = find_nearest(spectra_dict[key]['mu Ref'].energy, Ref_E0_calib-20)[0]
        ind_max = find_nearest(spectra_dict[key]['mu Ref'].energy, Ref_E0_calib+20)[0]
        larch.xafs.find_e0(spectra_dict[key]['mu Ref'].energy[ind_min:ind_max], 
                           mu = spectra_dict[key]['mu Ref'].mu[ind_min:ind_max], 
                           group = spectra_dict[key]['mu Ref'])
        Ref_E0_pos.append(spectra_dict[key]['mu Ref'].e0)
        
    
    print('Ref Calibraiton Statistics:')
    print('Ref E0 min: {0}'.format(min(Ref_E0_pos)))
    print('Ref E0 max: {0}'.format(max(Ref_E0_pos)))
    print('Ref E0 mean: {0}'.format(np.asarray(Ref_E0_pos).mean()))
    print('Ref E0 std: {0}'.format(np.asarray(Ref_E0_pos).std()))
    print('Ref E0 calibrated to: {0}\n\n'.format(Ref_E0_calib))


    #Shift energy scale based upon E0 values
    del_E = Ref_E0_calib - np.asarray(Ref_E0_pos).mean()

    for key in spectra_dict.keys():
        spectra_dict[key]['mu Ref'].energy = spectra_dict[key]['mu Ref'].energy + del_E
        spectra_dict[key]['mu Samp'].energy = spectra_dict[key]['mu Samp'].energy + del_E
    
    # Find new Edge Energy of each Reference Spectra & Sample Spectra


    Samp_E0_pos = []

    for key in spectra_dict.keys():
        Ref_ind_min = find_nearest(spectra_dict[key]['mu Ref'].energy, Ref_E0_calib-20)[0]
        Ref_ind_max = find_nearest(spectra_dict[key]['mu Ref'].energy, Ref_E0_calib+20)[0]
    
        Samp_ind_min = find_nearest(spectra_dict[key]['mu Samp'].energy, Samp_E0_calib-20)[0]
        Samp_ind_max = find_nearest(spectra_dict[key]['mu Samp'].energy, Samp_E0_calib+20)[0]
    
        larch.xafs.find_e0(spectra_dict[key]['mu Ref'].energy[Ref_ind_min:Ref_ind_max],
                       mu = spectra_dict[key]['mu Ref'].mu[Ref_ind_min:Ref_ind_max],
                       group = spectra_dict[key]['mu Ref'])
        larch.xafs.find_e0(spectra_dict[key]['mu Samp'].energy[Samp_ind_min:Samp_ind_max], 
                       mu = spectra_dict[key]['mu Samp'].mu[Samp_ind_min:Samp_ind_max],
                       group = spectra_dict[key]['mu Samp'])
    
        Samp_E0_pos.append(spectra_dict[key]['mu Samp'].e0)
    
    print('Sample Edge Energy Statistics:')
    print('Samp E0 min: {0}'.format(min(Samp_E0_pos)))
    print('Samp E0 max: {0}'.format(max(Samp_E0_pos)))
    print('Samp E0 mean: {0}'.format(np.asarray(Samp_E0_pos).mean()))
    print('Samp E0 std: {0}'.format(np.asarray(Samp_E0_pos).std()))
    
    return

def cluster2feffinp(atoms, absorbing_atom_index, edge, title, output_file_name, distance_cutoff = 8.0,is_DFT = False, DFT_rep = [2,2,1], r_max = 5.0):
    '''
    Takes an atoms cluster, often from a cif file, and saves an feff6 input file based upon a clsuter size wiht radial distance of "distance_cutoff"

    atoms : ase Atoms
        atoms object to generate cluster and feff input file on
    absorbing_atom_index : int
        index of atom in atoms to be used as the absorber and "center" of cluster
    edge : str
        X-ray edge that feff will run at, e.g. K, L1, L2...
    title : str
        what the title line will say in the feff inp file. could be ''
    output_file_name : str
        path to save file at
    distance_cutoff : float, optional
        longest range the cluster will simulate out to. passed to atoms2cluster 
        The default is 8.0.
    is_DFT : bool, optional
        True if the atoms object is to nto be scaled in one or more direction. passed to atoms2cluster. 
        The default is False.
    DFT_rep : list/array, optional
        limits that a defiend crystal is to be scaled. used in conjunction wiht is_def.
        passed to atoms2cluster. The default is [2,2,1].
    r_max : float, optional
        sets range to define feff6 calcualtion in the inp file. The default is 5.0.

    Returns
    -------
    cluster : ASE atoms object
        atoms object containing all atoms inside the cluster
        
    atoms_list : list
        list fo xyz coordiantes, ipot, tag, and radial distance from center atom.
        equivalent to the atoms table in the feff output file

    '''
    cluster, absorber_index = atoms2cluster(atoms, absorbing_atom_index, distance_cutoff, is_DFT = is_DFT , DFT_rep = DFT_rep)
    
    view(atoms)
    view(cluster)
    
    # Modified to be more generic
    title= title        
    core = cluster[absorber_index].symbol
    edge = edge
    HOLE = edge2hole(edge)
    edge_energy = xraylib.EdgeEnergy(xraylib.SymbolToAtomicNumber(core),HOLE-1)*1000
    #print(title, core, edge, HOLE, edge_energy)
    
    # Find the unique elements and their atomic numbers
    tag = np.unique(cluster.get_chemical_symbols())
    
    # Set teh first line of the potential list based upon the core    
    potentials = [[0, cluster.get_atomic_numbers()[absorber_index], cluster[absorber_index].symbol]]
    
    
    # Add lines to the potentail list based upon tags
    for i in range(len(tag)):
        potentials.append([i+1, xraylib.SymbolToAtomicNumber(tag[i]), tag[i]])
    print(potentials)
    
   
    atoms_list = []

    for atom in cluster:   
        temp_line = []
        
        temp_position = cluster.get_distance(absorber_index, atom.index, vector=True, mic=True).tolist()        
        temp_distance = cluster.get_distance(absorber_index, atom.index, mic= True)
        temp_symbol = cluster[atom.index].symbol
        
        if temp_distance == 0:
            temp_pot = 0
        else:
            for line in potentials[1:]:
                if temp_symbol == line[2]:
                    temp_pot = line[0]
                    
        temp_line.extend(temp_position)
        temp_line.append(temp_pot)
        temp_line.append(temp_symbol)
        temp_line.append(temp_distance)
        
        atoms_list.append(temp_line)
    
    atoms_list.sort(key = lambda x: x[-1])
    
    # Remove excess potentials from potentials list
    potentials2 = [potentials[0]]
    
    for i in range(1,len(potentials)):
        for line in cluster:
            if not line.index == absorber_index:
                if line.symbol == potentials[i][2] :
                    potentials2.append(potentials[i])
                    break
    
    
    file_text = f"""
 TITLE     {title}
 
 HOLE      {HOLE}   1.0  *  FYI: ({core} {edge} edge @ {edge_energy} eV, 2nd number is S0^2)
 *         mphase,mpath,mfeff,mchi
 CONTROL   1      1     1     1  
 PRINT     1      0     0     0  
 
 RMAX      {r_max:.1f}
  * POLARIZATION  0   0   0
 
 POTENTIALS
  * ipot   Z      tag
"""
    
    for line in potentials2:
        file_text += f"    {line[0]}       {line[1]}       {line[2]}\n"
        
    
    file_text += f"""
 ATOMS                  * this list contains {len(cluster)} atoms
 *  x              y              z        ipot      tag      distance
"""

    for line in atoms_list:
        file_text += f"{float(line[0]):11.5f}{float(line[1]):11.5f}{float(line[2]):11.5f}{int(line[3]):3d}       {line[4]}       {float(line[5]):.5f}\n"
    file_text += "END\n"
    with open (output_file_name, 'w') as file:
        file.write(file_text)

    print('Finished atoms2feff_inp')
    
    return cluster, atoms_list

def edge2hole(edge):
    '''
    A quick way to convert the string name of an edge to a hole value for feff6
    use in conjunction with cluster2feffinp

    Parameters
    ----------
    edge : str
        the edge that feff will calcualte scattering parameters at

    Returns
    -------
    hole : int
        interger value used by feff 6 to represent which excitation event.

    '''
    
    edges = ['K', 
             'L1', 'L2', 'L3', 
             'M1', 'M2', 'M3', 'M4', 'M5',
             'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7',
             'O1', 'O2', 'O3', 'O4', 'O5',
             'P1', 'P2', 'P3']
    
    hole = edges.index(edge) + 1
    
    return hole  

def feff_PathSummary(filename, sigma2 = 0.003):
    '''
    returns dictionary cotnaining feffpath groups and a summary of all paths 
    sorted by scattering type frm feff input file "filename"

    Parameters
    ----------
    filename : str
        full path to feff inpu file "\*.inp".
    sigma2 : float, optional
        dampening factor. The default is 0.003.

    Returns
    -------
    dictX : dictionary
        cotnaining feffpath groups and a summary of all paths sorted by
        back scatterer or multiple scattering "ms".

    '''
    
    # Creat subdirectory matching file name to store feff data, move feff inp file to directory
    fname = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    
    feff_dir = create_subdir(dirname, fname[:-4])
    feff_path = os.path.join(feff_dir,fname)
    
    copyfile(filename,feff_path)
    
      
    # Run feff6l on the feff input file 
    
    larch.xafs.feff6l(feff_path)
    
    # List all feff paths that were created
    feff_paths = glob2.glob(os.path.join(os.path.dirname(feff_path),'feff*.dat'))
    
    # Build path resutls dictionary
    dictX = {}
    
    # Get all paths into result dictionary
    for line in feff_paths:
        gp = larch.xafs.feffpath(line)
        larch.xafs.path2chi(gp)
        gp.sigma2 = sigma2
    
        if gp.nleg > 2:
            gp.label = 'ms'
        else:
            gp.label = gp.geom[1][0]
    
        dictX[os.path.basename(line)] = gp
    
    # Find Paths and extract key infomraiton
    
    path_bs = []
    
    path_reff = []
    
    path_degen = []
    
    for key in dictX.keys():
        path_bs.append(dictX[key].label)
        path_reff.append(dictX[key].reff)
        path_degen.append(dictX[key].degen)
    
    path_summary_df = pd.DataFrame({"feff_path":line[:-4], "Backsactterer": path_bs, "Reff": path_reff, "Degen": path_degen})
    
    # Split paths into dictionary for easy plotting
    unique_paths = path_summary_df.Backsactterer.unique()
    
    dict1 = {}
    
    for name in unique_paths:
        dict1[name] = path_summary_df[path_summary_df['Backsactterer'] == name]
    
    dictX["Path_Summary"] = dict1
        

    return dictX


def Interp_Spectra(spectra_dict, offset = 5):
    '''
    NEEDS UPDATING - based on dictionary from ReadCXAS

    Parameters
    ----------
    spectra_dict : TYPE
        DESCRIPTION.
    offset : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    None.

    '''
    # Empty lists to store energy values
    low_limit = []
    high_limit = []

    # Find the highest low and lowest high for all energy scales in the dictionary
    for key in spectra_dict.keys():
        low_limit.append(spectra_dict[key]['mu Samp'].energy.min())
        high_limit.append(spectra_dict[key]['mu Samp'].energy.max())
    
    print('Starting Energy Statistics:')
    print(f'\t Minimim Starting Energy: {min(low_limit)}')
    print(f'\t Maximum Starting Energy: {max(low_limit)}')
    print(f'\t Starting Energy mean: {np.asarray(low_limit).mean()}')
    print(f'\t Starting Energy std: {np.asarray(low_limit).std()}')
    
    print('\nEnding Energy Statistics:')
    print(f'\t Minimim Ending Energy: {min(high_limit)}')
    print(f'\t Maximum Ending Energy: {max(high_limit)}')
    print(f'\t Ending Energy mean: {np.asarray(high_limit).mean()}')
    print(f'\t Ending Energy std: {np.asarray(high_limit).std()}')
    
    # Find index of first spectra that are bounded by the offset
    offset = offset
    temp_key = next(iter(spectra_dict)) 
    min_ind = find_nearest(spectra_dict[temp_key]['mu Samp'].energy, max(low_limit)+offset)[0]
    max_ind = find_nearest(spectra_dict[temp_key]['mu Samp'].energy, min(high_limit)-offset)[0]
    
    print('\nValues of bound energies: {0}, {1}'.format(spectra_dict[temp_key]['mu Samp'].energy[min_ind],spectra_dict[temp_key]['mu Samp'].energy[max_ind]))
    

    # Determine the new index to interpoalte onto from the index limits above
    CommonEnergy_index = spectra_dict[temp_key]['mu Samp'].energy[min_ind:max_ind]

    #Interpolate all data sets onto new axis, and build an interpoalted dictionary
    for key in spectra_dict.keys():
        d = {'Sample': spectra_dict[key]['mu Samp'].mu,
             'Reference': spectra_dict[key]['mu Ref'].mu}
        index = spectra_dict[key]['mu Samp'].energy

        temp_df = pd.DataFrame(data = d, index = index)
        temp_df.index.rename('Energy', inplace = True)
        
        interpEnergy_df = interp_df(temp_df, CommonEnergy_index)
        
        spectra_dict[key]['mu Samp'].energy = CommonEnergy_index
        spectra_dict[key]['mu Samp'].mu = interpEnergy_df['Sample'].values
        spectra_dict[key]['mu Ref'].energy = CommonEnergy_index
        spectra_dict[key]['mu Ref'].mu = interpEnergy_df['Reference'].values
        
        d.clear()
    return

def Normalize_Spectra(spectra_dict, pre1, pre2, norm1, norm2, nnorm, make_flat):
    '''
    FILL ME IN

    Parameters
    ----------
    spectra_dict : TYPE
        DESCRIPTION.
    pre1 : TYPE
        DESCRIPTION.
    pre2 : TYPE
        DESCRIPTION.
    norm1 : TYPE
        DESCRIPTION.
    norm2 : TYPE
        DESCRIPTION.
    nnorm : TYPE
        DESCRIPTION.
    make_flat : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    for key in spectra_dict.keys():
        larch.xafs.pre_edge(spectra_dict[key]['mu Ref'].energy, mu = spectra_dict[key]['mu Ref'].mu, 
                            group = spectra_dict[key]['mu Ref'], e0 = spectra_dict[key]['mu Ref'].e0, 
                            pre1 = pre1, pre2 = pre2, norm1 = norm1, norm2 = norm2, 
                            nnorm = nnorm, make_flat = make_flat)
    
        larch.xafs.pre_edge(spectra_dict[key]['mu Samp'].energy, mu = spectra_dict[key]['mu Samp'].mu, 
                            group = spectra_dict[key]['mu Samp'], e0 = spectra_dict[key]['mu Samp'].e0,
                            pre1 = pre1, pre2 = pre2,  norm1 = norm1, norm2 = norm2, 
                            nnorm = nnorm, make_flat = make_flat)
    
    return

def plot_feff_paths(feff_dict, xmin= 0, xmax=6, xlabel = 'R (Å)', ylabel = 'Coordination Number'):
    '''
    quick visualizatiaon of scattering path distributions from a feff path summary
    used in conjunction wiht dictionary outout from feff_PathSummary
    
    figure will be displayed in line.

    Parameters
    ----------
    feff_dict : dictionary
        dectionary generated from feff_PathSummary.
    xmin : int/float, optional
        x axis lower limit. The default is 0.
    xmax : int/float, optional
        x axis upper limit. The default is 6.
    xlabel : str, optional
        y axis label. The default is 'R (Å)'.
    ylabel : str, optional
        y axis labe. The default is 'Coordination Number'.

    Returns
    -------
    matplotlib figure

    '''    
    fig1 = plt.figure(constrained_layout=True, figsize = (10,10))
    spec1 = gridspec.GridSpec(ncols = 1, nrows = 1, figure = fig1)

    f1_ax1 = fig1.add_subplot(spec1[0,0])
    
    path_dict = feff_dict['Path_Summary']
    
    for key in path_dict:
        f1_ax1.scatter(path_dict[key].Reff.values, path_dict[key].Degen.values, label = key)
        f1_ax1.vlines(path_dict[key].Reff.values, 0, path_dict[key].Degen.values, color = 'k')
    
    f1_ax1.legend(loc = 'upper left')
    f1_ax1.set_xlim(xmin,xmax)
    f1_ax1.set_xlabel(xlabel)
    f1_ax1.set_ylabel(ylabel)
    
    return fig1
    
def ReadCXAS(files, Energy_col, muSamp_N, muSamp_D, muRef_N, muRef_D, 
             muSamp_state = True, muSamp_flip = False, 
             muRef_state = True, muRef_flip = False, 
             time_stamp = False, padded = False):
    """
    Reads in a list of CXAS data files (.txt), calculates the absorption 
    coefficient of the sample and the refernece, and stores each as a larch
    group in a dictionary keyed to each samples file name.
    
    **Note** Current files must be from the database CXAS data collector.
    Other file formats are not currently supported

    Parameters
    ----------
    files : list
        List of files (wiht full path) to read.
    Energy_col : int
        column in the data file corresponding to enegy.
    muSamp_N : int
        column in the data file corresponding to the numerator 
        used to calculate the absorption coefficient of the sample.
    muSamp_D : int
        column in the data file corresponding to the denominator 
        used to calculate the absorption coefficient of the sample.
    muRef_N : int
        column in the data file corresponding to the numerator 
        used to calculate the absorption coefficient of the reference.
    muRef_D : int
        column in the data file corresponding to the denominator 
        used to calculate the absorption coefficient of the reference.
    muSamp_state : bool, optional
        used to set log state in calc_mu() for sample. If True, ln() applied.
        The default is True.
    muSamp_flip : bool, optional
        used to set flip state in calc_mu() for sample. If False, N/D applied.
        The default is False.
    muRef_state : bool, optional
        used to set log state in calc_mu() for reference. If True, ln() applied.
        The default is True.
    muRef_flip : bool, optional
        used to set flip state in calc_mu() for reference. If False, N/D applied.
        The default is False.
    time_stamp : bool, optional
        used to set time or spectra index in CXAS_Sorted(). If True, index is 
        sorted by time spectra was collected. The default is False.
    padded : bool, optional
        used to define if spectra numbers are padded in CXAS_Sorted(). If False,
        spectra file names are padded on the order of XXXX. The default is False.

    Returns
    -------
    Spectra : Dictionary
        Dictionary with each key as a spectra file name. Each key contains
        two populated larch groups, one for mu Sample, the other for mu Reference.
    file_list : DataFrame
         equivalent DataFrame to thot of CXAS_Sorted(). Indexed to either time
         or scan number, containing file names, padded names, and file paths

    """
    
    # Sort, and possibly pad, filenames by colelction time or spectra number
    file_list = CXAS_Sorted(files, time_stamp = time_stamp, padded = padded)
    
    # Result Dictionary
    Spectra = {}

    # Populate dictionary
    for index, row in file_list.iterrows():
        filename = row['Padded Name']
    
        #Build Dictionary containing dictionary of file names, with Larch goups for the sample and the reference    
        Spectra[filename] = {}
        Spectra[filename]['Time'] = index 
        Spectra[filename]['mu Samp'] = larch.Group(name = str(filename + '_mu Samp'))
        Spectra[filename]['mu Ref'] = larch.Group(name = str(filename + '_mu Ref'))

        # Extract Detector Channles and store in array:
        raw_data = np.genfromtxt(row['Path'])
    
    
        # Make df of signals and energy, then remove any NaN values
        energy = raw_data[:, Energy_col]
        mu_samp = calc_mu(raw_data[:, muSamp_N], raw_data[:, muSamp_D], muSamp_state, flip = muSamp_flip)
        mu_ref = calc_mu(raw_data[:, muRef_N], raw_data[:, muRef_D], muRef_state, flip = muRef_flip)
    
        # Create a temporary dictionary of energy and adsorption coefficients for ease of making DataFrame    
        temp_data = {'Energy': energy, 'mu Samp': mu_samp, 'mu Ref': mu_ref}
    
        #Remove +/-inf and NaN rows
        temp_df = pd.DataFrame(temp_data)
        temp_df.sort_values(by = ['Energy'], inplace = True)
        temp_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        temp_df.dropna(axis=0,inplace = True)
   
    
        # Add Energy to Larch Groups:    
        Spectra[filename]['mu Samp'].energy = temp_df['Energy'].values
        Spectra[filename]['mu Ref'].energy = temp_df['Energy'].values
    
        # Add mu Sample to larch Group:
        Spectra[filename]['mu Samp'].mu = temp_df['mu Samp'].values
    
        # Add mu Ref to Larch Group:
        Spectra[filename]['mu Ref'].mu = temp_df['mu Ref'].values
    
    return Spectra, file_list


def ReadFTIR(directory):
    """
    Combines a directory of csv FTIR spectra with timestamps included from the RUN SAMPLE macro.
    
    Retuns a pandas datframe wiht index = wavenumber, and columns are spetra wiht timestampe names.
    """
    
    #List all the .csv files in the folder
    files = glob2.glob(directory + '/*.csv')
    
    #Time format in file name:
    fmt = "%b %d %H-%M-%S %Y"
    
    #Determine nuber of files
    data_files = np.shape(files)[0]
    
    # Extract wave number and adsorption spectra for each file
    for i in range(data_files):

        # Store timestamp string at end of file
        time = os.path.basename(files[i])[-24:-4]
        
        if i == 0:
            #Convert to dataframe
            data_df = pd.read_csv(files[i], sep = ',', names = ['wavenumber', dt.strptime(time,fmt)])
            data_df.set_index('wavenumber', inplace = True)
        else:
            #Convert to dataframe
            data_tmp = pd.read_csv(files[i], sep = ',', names = ['wavenumber', dt.strptime(time,fmt)])
            data_tmp.set_index('wavenumber', inplace = True)
            
            data_df = pd.concat([data_df, data_tmp], axis = 1)
        
    return data_df


def Set_t_Zero(date_str, date_fmt, df):
    """
    date_str: when time = 0 in experiment.
    date_fmt = time format; suggest '%m/%d/%Y %I:%M:%S %p'
    df: pandas dataframe to shift time index
    
    retuns df with new time index
    """
    temp_df = df
    
    dt_time = dt.strptime(date_str, date_fmt)
    
    times = temp_df.index
    #print(times)
    
    minutes = (times - dt_time).total_seconds() // 60
    #print(minutes)
    
    temp_df['tZero'] = minutes
    
    temp_df.reset_index(inplace = True)
    #print(temp_df)
    
    temp_df.set_index('tZero', inplace = True)
    #print(temp_df)
    
    return temp_df