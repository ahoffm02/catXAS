# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:14:21 2022

@author: ashoff
"""

##############################################################################

                                # Modules #
                        
##############################################################################

# File Handling
import os
import glob2 as glob
import pickle

# Timestamps
import datetime
from datetime import datetime as dt

# Data organization
import pandas as pd
import numpy as np

# X-ray Science
import larch
from larch.io import read_ascii
from larch.io import merge_groups


# Plotting
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

# From Catxas
import general as fcts
import xas as xfcts
import plot as pfcts
import process


##############################################################################

            # EXPERIMENTAL CLASS Operating Functions #
                        
##############################################################################

def open_experiment(fname):
    '''
    Opens and experimental class object using pickling

    Parameters
    ----------
    fname : STR
        full path and file name, extension optional to the .pickle file.

    Returns
    -------
    my_exp : experiment object
        experimental class object.

    '''
    if fname[-7:] == '.pickle':
        pass
    else:
        fname = fname+'.pickle'
        
    pickle_in = open(fname,"rb")
    my_exp = pickle.load(pickle_in)
    pickle_in.close()

    return my_exp

def save_experiment(experiment, fname):
    '''
    Saves and experimental class object using pickling

    Parameters
    ----------
    experiment : experiment object
        experiment class object.
    fname : STR
        Path and filename where to save the experimetnal object. Extension optional.

    Returns
    -------
    None.

    '''
    if fname[-7:] == '.pickle':
        pass
    else:
        fname = fname+'.pickle'

    pickle_out = open(fname,"wb")
    pickle.dump(experiment, pickle_out)
    pickle_out.close()
    
    return

def groups_lists(exp, spectra_name_list, spectra_name = 'mu sample'):
    grp_list = []
    for key in exp.spectra.keys():
        if key in spectra_name_list:
            grp_list.append(exp.spectra[key]['Absorption Spectra'][spectra_name])
    
    return grp_list


def process_concat(exp, spectra_name_list):
    process_df_list = []
    for key in exp.spectra.keys():
        if key in spectra_name_list:
            process_df_list.append(exp.spectra[key]['Process Values'])
            
    process_df = pd.concat(process_df_list)
    
    process_df.drop(columns=['File Name'], inplace = True)
    
    process_df.reset_index(drop = True, inplace = True)
    
    return process_df

def time_lists(exp, spectra_name_list):
    time_list = []
    for key in exp.spectra.keys():
        if key in spectra_name_list:
            time_list.append(exp.spectra[key]['Time'])
            
    return time_list

def merge_spectra(exp, spectra_name_list, xarray = 'energy', yarray = 'mu'):
    
    # Create List of Groups [Sample and Ref]
    S_grp_list = groups_lists(exp, spectra_name_list, spectra_name = 'mu Sample')
    R_grp_list = groups_lists(exp, spectra_name_list, spectra_name = 'mu Reference')
    
    
    # Create list of Process Parameters
    process_params = process_concat(exp, spectra_name_list)
    
    # Create list of Timestamps
    time_list = time_lists(exp, spectra_name_list)
    
    # Calculate the mean absorption spectra [Sample and Ref]
    S_grp_mean = merge_groups(S_grp_list, master = None, 
                          xarray = xarray, yarray = yarray, 
                          kind = 'cubic', trim = True, calc_yerr = True)
    
    R_grp_mean = merge_groups(R_grp_list, master = None, 
                          xarray = xarray, yarray = yarray, 
                          kind = 'cubic', trim = True, calc_yerr = True)
    
    # Calculate the mean Process Paramters
    process_mean = process_params.mean(axis = 0)
    
    # Calculate the mean Timestamp
    time_mean = pd.Timestamp(np.mean([i.timestamp() for i in time_list]), unit = 's')
    
    # Rename group based upon TOS [Sample and Ref]
    S_grp_mean.__name__ = f'{process_mean["TOS [s]"]:0.2f} s'
    R_grp_mean.__name__ = f'{process_mean["TOS [s]"]:0.2f} s - Ref'
    
    
    # Create Spectra dictionary for experiment class
    spectra_dict = {'XAS Data Structure': None, 
                    'Time': time_mean, 
                    'Absorption Spectra':{'mu Sample': S_grp_mean, 'mu Reference': R_grp_mean},
                    'Process Values':process_mean
                   }
    
    return spectra_dict

##############################################################################

            # EXPERIMENTAL CLASS - IN-SITU XAS PROCESSING #
                        
##############################################################################

class Experiment:
    
    ########################
    
    ###### Constructor #####
    
    ########################
    
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
         # Update 1/28/2024 by ASH - Force empty XAS Spectra Files datafame
        self.summary = {'XAS Spectra Files': pd.DataFrame(columns=['Time','TOS [s]', 'File Name', 'Padded Name', 'Path'])}
        self.summary['XAS Spectra Files'].set_index('Time', inplace=True)
        
        return
    
    ##########################################################################
    
    ############################ Import Functions ############################
    
    ##########################################################################
    
    ##############################################
    
    ##### Import Process Parameter Functions #####
    
    ##############################################
    
    def import_massspec(self, file_name):
        '''
        Loads mass spectrometer data from Hiden Mass Specs. into the 
        process_params dictionary of the experimental class

        Parameters
        ----------
        file_name : STR
            Full path string to mass spec. file.

        Returns
        -------
        None.

        '''
        
        self.process_params['MS Data'] = process.ReadMSData(file_name)
        self.summary['MS Filename'] = file_name
        
        return
        
    
    def import_labview(self, file_name):
        '''
        Loads LabView data generated from the Co-ACCESS flow systems into the 
        process_params dictionary of the experimental class

        Parameters
        ----------
        file_name : STR
            Full path string to LabView file.

        Returns
        -------
        None.

        '''
        
        self.process_params['LV Data'] = process.ReadLVData(file_name)
        self.summary['LV Filename'] = file_name
        
        return

    
    def import_process_df(self, process_df, name = 'process_data'):
        '''
        Uploads process data in the form of a time-indexed datafram

        Parameters
        ----------
        process_df : Pandas Series or Dataframe with a datetime index
            Pandas Dataframe or Series with a datetime index. Columns represent
            process data recorded at each time stamp.
        name : str, optional
            A name to represent the type of data stored in the dataframe. The default is 'process_data'.

        Returns
        -------
        None.

        '''
        
        # Check if name is already used
        if name in self.process_params.keys():
            
            print('Name of Process data alerady exists. Choose a new name.')
            
        # Check if df index is in a datetime format
        elif isinstance(process_df.index, pd.DatetimeIndex) == False:
            
            print('Index of data is not in a datetime format. Redefine the inext and upload again.')
            
        # Store the data
        else:
            self.process_params[name] = process_df
            self.summary[name] = f'"{name}" data uploaded as a dataframe and not from a file'
            
        return
        
    ##############################################
    
    #####      Import Spectra Functions      #####
    
    ##############################################
    
    def import_spectrum(self, file, xas_data_structure, print_name = False):
        #Added 1/27/2024 by ASH - goal is to creaate method to append filed to an empty or existing experiment
        '''
        Loads XAS data of a single spectrum into the spectra dictionary of the experiment class.
        
        Creates a time/file/path dataframe in the 'summary' dictionary of the
        experiment class.

        Parameters
        ----------
        fname : STR
            Path to single spectrum file.
        xas_data_structure : DICT
            Dictionary of XAS data strcuture. Examples found in 
            "BL specific XAS data structures.ipynb" in teh notebooks section
            of the module.
        print_name : BOOL, optional
            Returns the name of each spectra uploaded to experiment object.
            Used for troubleshooting. The default is False.

        Returns
        -------
        None.

        '''        
        
        # Get parameters from the data structure
        
        time_stamp = xas_data_structure['time stamp']
        time_line = xas_data_structure['time on line']
        time_format = xas_data_structure['time format']
        padded = xas_data_structure['padded scan numbers']
        is_QEXAFS = xas_data_structure['is QEXAFS']
        col_names = xas_data_structure['column names']
        energy_name = xas_data_structure['energy column']
        
        # Get the name of the scan without path or extension
        fname = os.path.basename(file)[:-4]
        
        # Added padded zeros to scan if sample does not have any
        if not padded:
            scan_no = fcts.get_trailing_number(fname)
            char = len(str(scan_no))
            fname_padded = fname[:-char] + str(scan_no).zfill(4)
        else:
            fname_padded = fname
        
        # Display name if the flag is called
        if print_name:
            print(fname)

        # Build the scan dictionary in the spectra dictionary
        self.spectra[fname] = {}
    
        self.spectra[fname]['XAS Data Structure'] = xas_data_structure
        
        self.spectra[fname]['BL Data'] = read_ascii(file, labels=col_names)
        
        self.spectra[fname]['BL Data'].__name__ = fname
        
        # Determine the enregy values at the start/end of the dataset
        Einit = self.spectra[fname]['BL Data'].__dict__[energy_name][0]
        Efin = self.spectra[fname]['BL Data'].__dict__[energy_name][-1]
        
        # Invert the dataset if the energy values are not increasing
        if Einit >= Efin:
            for line in col_names:
                self.spectra[fname]['BL Data'].__dict__[line] = np.flipud(self.spectra[fname]['BL Data'].__dict__[line])

        
        
        # Updating the Summary File
        if time_stamp:
            
            time_str = self.spectra[fname]['BL Data'].header[time_line]
            
            if not is_QEXAFS:
                time = dt.strptime(time_str, time_format)

            elif is_QEXAFS:
                # To be checked with new QXAFS Data
                time = dt.strptime(time_str, time_format)
                micros = self.spectra[fname]['BL Data'].time[0] # to be checked
                micros_to_date = datetime.timedelta(microseconds = micros)
                time = time + micros_to_date

            else:
                pass
       
            temp_dict = {'Time': [time],
                         'TOS [s]': [0],
                         'File Name': [fname],
                         'Padded Name': [fname_padded],
                         'Path': [file]}
            
            temp_df = pd.DataFrame(temp_dict)
            temp_df.set_index('Time', inplace = True)
            
        else:
            temp_dict = {'Time': [dt.time()],
                         'TOS [s]': [0],
                         'File Name': [fname],
                         'Padded Name': [fname_padded],
                         'Path': [file]}
            
            temp_df = pd.DataFrame(temp_dict)
            temp_df.set_index('Time', inplace = True)
                        
        self.summary['XAS Spectra Files'] = pd.concat([self.summary['XAS Spectra Files'], temp_df], axis = 0, ignore_index=False)
                
        # Add Time to Sample
        self.spectra[fname]['Time'] = time
        
        # Rearrange Summary File and Spectra Dictionary. Add TOS if appropriate
        if time_stamp:
            self.summary['XAS Spectra Files'].sort_index(axis=0, inplace=True)
            self.spectra = {i: self.spectra[i] for i in self.summary['XAS Spectra Files']['File Name'].values}
            
            # Determine the time elapses since the first file in seconds
            elapsed_time = []
            for line in self.summary['XAS Spectra Files'].index.values:
                temp_TOS = line-self.summary['XAS Spectra Files'].index.values[0]
                temp_TOS = temp_TOS / np.timedelta64(1, 's')
                elapsed_time.append(temp_TOS)
            self.summary['XAS Spectra Files']['TOS [s]'] = elapsed_time
            
            
        else:
            self.summary['XAS Spectra Files'].sort_values(by = 'Padded Name', ignore_index = True, inplace = True)
            self.spectra = {i: self.spectra[i] for i in self.summary['XAS Spectra Files']['File Name'].values}
            
        
        return
    
    def import_spectra(self, xas_data_directory, xas_data_structure, ext = '.txt', print_name = False):
        '''
        Function to import multiple spectra from a directory (or subfolders in a directory) using the simplified import_spectrum fucntion

        Parameters
        ----------
        xas_data_directory : STR
            Path to XAS spectra to be imported.
        xas_data_structure : DICT
            Dictionary of data structure to be loaded into the experiment with each file.
        ext : STR, optional
            extension of the spectra files for importing. The default is '.txt'.
        print_name : TYPE, optional
            Prints the spectra name after upload to help with debugging. The default is False.

        Returns
        -------
        None.

        '''
        
        # Use glob2 to get a list of all files in files_directory
        files = glob.glob(xas_data_directory+f'/**/*{ext}', recursive=True)
        
        for line in files:
        
            self.import_spectrum(line, xas_data_structure, print_name = print_name)
        
        return
    
    
    def import_spectra_data(self, xas_data_directory, xas_data_structure, print_name = False):
        '''
        Depreciated. See import_spectra or import_spectrum.
        
        Loads XAS data into the spectra dictionary of the experiment class.
        
        Creates a time/file/path dataframe in the 'summary' dictionary of the
        experiment class.

        Parameters
        ----------
        xas_data_directory : STR
            Path to diretory where XAS data is stored.
        xas_data_structure : DICT
            Dictionary of XAS data strcuture. Examples found in 
            "BL specific XAS data structures.ipynb" in teh notebooks section
            of the module.
        print_name : BOOL, optional
            Returns the name of each spectra uploaded to experiment object.
            Used for troubleshooting. The default is False.

        Returns
        -------
        None.

        '''  
        
        print('Depreciated, use import_spectra for multiple files, or import_spectrum for a single file')
        
        return
    
       
    #################################
    
    ##### Correlation Functions #####
    
    #################################
    
    def correlate_process_params(self):
        '''
        TBD
        '''          
        
        self.summary['XAS Spectra Process Params'] = self.summary['XAS Spectra Files'][['File Name','TOS [s]']]
        print('Genearted Spectra Summary')
        
        for key in self.process_params.keys():                
            temp_df = fcts.mergeindex(self.summary['XAS Spectra Process Params'], self.process_params[key])
            print(f'Merged Index for Process: {key}')
            
            self.summary['XAS Spectra Process Params'] = pd.concat([self.summary['XAS Spectra Process Params'], temp_df], axis=1)
            
        for key in self.spectra.keys():
            self.spectra[key]['Process Values'] =  self.summary['XAS Spectra Process Params'].loc[[self.spectra[key]['Time']]]
        
        return
    
    ###########################################
    
    ##### Spectra Interrogation Functions #####
    
    ###########################################
    
    def data_length_screen(self, deviations = 3, print_summary = True):
        # Extract the number of data points and the starting and ending E value in each spectrum for interrogation
        fnames = []
        datapts = []
        Estart = []
        Eend = []
                
        for key in self.spectra.keys():
            fnames.append(key)
            datapts.append(len(self.spectra[key]['BL Data'].Energy))
            Estart.append(self.spectra[key]['BL Data'].Energy[0])
            Eend.append(self.spectra[key]['BL Data'].Energy[-1])
            
        data_dict = {'Filename': fnames,
                   'Data Points': datapts,
                   'Start Energy': Estart,
                   'End Energy': Eend}
        
        data_df = pd.DataFrame(data_dict) 

        # Determine the statistical characteristics of the data points
        mean_dpts = data_df['Data Points'].mean()
        stdev_dpts = data_df['Data Points'].std()
        
        # Determine the statistical characteristics of the starting energy values
        mean_Estart = data_df['Start Energy'].mean()
        stdev_Estart = data_df['Start Energy'].std()
        
        # Determine the statistical characteristics of the starting energy values
        mean_Eend = data_df['End Energy'].mean()
        stdev_Eend = data_df['End Energy'].std()
        
        # Set Threshold limits
        low_pts = mean_dpts-deviations*stdev_dpts
        high_pts = mean_dpts+deviations*stdev_dpts
        
        low_Estart = mean_Estart-deviations*stdev_Estart
        high_Estart = mean_Estart+deviations*stdev_Estart
        
        low_Eend = mean_Eend-deviations*stdev_Eend
        high_Eend = mean_Eend+deviations*stdev_Eend
        
        # extract all spectra that are outside of the data point threshold 
        problem_spectra = data_df[(data_df['Data Points']<low_pts) | 
                                     (data_df['Data Points']>high_pts) | 
                                     (data_df['Start Energy']<low_Estart) | 
                                     (data_df['Start Energy']>high_Estart) | 
                                     (data_df['End Energy']<low_Eend) | 
                                     (data_df['End Energy']>high_Eend)]
        
        # Write out what was tested and what as found
        if print_summary:
            print("\u0332".join("Spectra Data Length/Starting Energy Characteristics:"))
            print('\n')
            print(f'\tSpectra interrogated: {len(data_df)}')
            print('\n')
            print(f'\tLongest Data Set: {data_df["Data Points"].max()} data points')
            print(f'\tShortest Data Set: {data_df["Data Points"].min()} data points')
            print(f'\tMean Data Points per Spectrum: {mean_dpts:0.0f}')
            print(f'\tDeviation in Data Points: {stdev_dpts:0.0f}')
            print('\n')
            print(f'\tLargest Starting Energy: {data_df["Start Energy"].max()} eV')
            print(f'\tSmallest Starting Energy: {data_df["Start Energy"].min()} eV')
            print(f'\tMean Starting Enregy: {mean_Estart:0.0f}')
            print(f'\tDeviation in Starting Energy: {stdev_Estart:0.0f}')
            print('\n')
            print(f'\tLargest Starting Energy: {data_df["End Energy"].max()} eV')
            print(f'\tSmallest Starting Energy: {data_df["End Energy"].min()} eV')
            print(f'\tMean Starting Enregy: {mean_Eend:0.0f}')
            print(f'\tDeviation in Starting Energy: {stdev_Eend:0.0f}')
            print('\n')
            
            
            print("\u0332".join("Problematic Spectra:"))
            print(f'\tNumber of Spectra: {len(problem_spectra)}')
            for index, row in problem_spectra.iterrows():
                for index, value in row.items():
                    if index == 'Filename':
                        print(f'\t{(index)}: {value}')
                    else:
                        if type(value) == str or type(value) == int:
                            print(f'\t\t{(index+":").rjust(20)}\t{value}')
                        elif type(value) == float:
                            print(f'\t\t{(index+":").rjust(20)}\t{value:0.2f}')
            
        # Returns df of bad spectra with their edge steps
        return problem_spectra  
            
    
    def edge_step_screen(self, deviations = 3, print_summary = True, show_problem_spectra = True):
        # Extract all the edge steps for interrogation
        fnames = []
        steps = []
    
        for key in self.spectra.keys():
            fnames.append(key)
            steps.append(self.spectra[key]['Absorption Spectra']['mu Sample'].edge_step)
    
        edge_step_df = pd.concat([pd.Series(fnames, name='Filename'), pd.Series(steps, name='Edge Step')], axis = 1)
        
        # Determine the statistical characteristics of the edge steps
        mean_step = edge_step_df['Edge Step'].mean()
        stdev_step = edge_step_df['Edge Step'].std()
        
        # Set Threshold limits
        low_step = mean_step-deviations*stdev_step
        high_step = mean_step+deviations*stdev_step
        
        # extract all spectra that are outside of the edge step threshold 
        problem_spectra = edge_step_df[(edge_step_df['Edge Step']<low_step) | (edge_step_df['Edge Step']>high_step)]
        
        # Write out what was tested and what as found
        if print_summary:
            print("\u0332".join("Edge Step Characteristics:"))
            print(f'\tSpectra interrogated: {len(edge_step_df)}')
            print(f'\tLargest Edge Step: {edge_step_df["Edge Step"].max():0.3f}')
            print(f'\tSmallest Edge Step: {edge_step_df["Edge Step"].min():0.3f}')
            print(f'\tMean Edge Step: {mean_step:0.3f}')
            print(f'\tDeviation in Edge Step: {stdev_step:0.3f}')
            print('\n')
            print("\u0332".join("Problematic Spectra:"))
            print(f'\tNumber of Spectra: {len(problem_spectra)}')
            for index, row in problem_spectra.iterrows():
                    print(f'\t\t{row.to_string(header=False, index=False)}')
            
        # Plot bad spectra if requested
        if show_problem_spectra:
            for line in problem_spectra.Filename:
                x = self.spectra[line]['Absorption Spectra']['mu Sample'].energy
                y = self.spectra[line]['Absorption Spectra']['mu Sample'].mu
                plt.plot(x,y,label = line)
                if len(problem_spectra) <=10:                
                    plt.legend()
            
        # Returns df of bad spectra with their edge steps
        return problem_spectra, edge_step_df # edge_step_df is new, see if this breaks somewhere else (ASH 1/29/2024)
    
    def remove_bad_spectra(self, spectra_list, test_removal = True):
        '''
        Removes spectra files contained in spectra_list from (1) self.spectra, (2) lines form self.summary['XAS Spectra Files'], 
        and (3) self.summary['XAS Spectra Process Params'] (if it exists)
        
        Adds/appends a list in self.summary to indicate spectra removed in this step.
        
        Parameters
        ----------
        spectra_list : LIST
            List of filenames to remove from experiment class

        Returns
        -------
        None.

        '''
        # Test to see if Process Params Summary df exists
        try:
            self.summary['XAS Spectra Process Params']
        except KeyError:
            ParamSummary_exists = False
        else:
            ParamSummary_exists = True
        
    
        for line in spectra_list:
            print(f'Removing {line}')
            
            # Removes the spectra from the *.spectra dictionary
            self.spectra.pop(line, None)
            
            # Removes the spectra from the 'XAS Spectra Files' dataframe
            self.summary['XAS Spectra Files'] = self.summary['XAS Spectra Files'][self.summary['XAS Spectra Files']['File Name'] != line]
            
            # Removed the spectra from the 'XAS Spectra Process Params' dataframe if it exists
            if ParamSummary_exists:
                self.summary['XAS Spectra Process Params'] = self.summary['XAS Spectra Process Params'][self.summary['XAS Spectra Process Params']['File Name'] != line]
            
            # Test to see if the spectrum has been removed from the spectra list and 1-2 summary files:
            if test_removal:
                try:
                    self.spectra[line]
                except KeyError:
                    print('\tRemoved from Spectra')
                else:
                    pass
        
                try:
                    self.summary['XAS Spectra Files'][line]
                except KeyError:
                    print('\tRemoved from the summary: "XAS Spectra Files')
                else:
                    pass
        
                if ParamSummary_exists:         
                    try:
                        self.summary['XAS Spectra Process Params'][line]
                    except KeyError:
                        print('\tRemoved from the summary: "XAS Spectra Process Params"')
                    else:
                        pass   

        if 'Spectra Removed' in self.summary:
            self.summary['Spectra Removed'].append(spectra_list)
            print('Spectra Removed list has been updated')
        else:
            self.summary['Spectra Removed'] = spectra_list
            print('Spectra Removed list has been created in Summary')
    
    
    def check_Energy_Range(self, spectra_name = 'mu Sample', has_e0 = False, print_summary = True):
        
        # Build dataframe indlucing e0
        if has_e0:
            column_names = ["E_min", "E_max", "E0", "Min_E_Step", "Max_E_Step", "Mean_E_Step", 'STD_E_Step']
            df = pd.DataFrame(columns = column_names)
        
        # Build dataframe without 30
        if not has_e0:
            column_names = ["E_min", "E_max", "Min_E_Step", "Max_E_Step", "Mean_E_Step", 'STD_E_Step']
            df = pd.DataFrame(columns = column_names)
    
        # Populate dataframe with values from each spectrum    
        for key in self.spectra.keys():
            
            # Min and max energy of the spectra
            emin = min(self.spectra[key]['Absorption Spectra'][spectra_name].energy)
            emax = max(self.spectra[key]['Absorption Spectra'][spectra_name].energy)
            
            # Attempt to find e0
            try:
                e0 = self.spectra[key]['Absorption Spectra'][spectra_name].e0
            except Exception:
                e0 = None
            
            # Find the energy step between each data point
            temp_step = []
            for i in range(len(self.spectra[key]['Absorption Spectra'][spectra_name].energy)-1):
                temp_step.append(self.spectra[key]['Absorption Spectra'][spectra_name].energy[i+1]-self.spectra[key]['Absorption Spectra'][spectra_name].energy[i])
            
            # Calculate energy step based parameters 
            min_step = min(temp_step)
            max_step = max(temp_step)
            mean_step = np.asarray(temp_step).mean()
            std_step= np.asarray(temp_step).std()  
        
            # Populate Databases
            if has_e0:
                df2 = pd.DataFrame([[emin, emax, e0, min_step, max_step, mean_step, std_step]], columns=column_names)
                
                df = pd.concat([df, df2], axis = 0, ignore_index = True)
                
            if not has_e0:
                df2 = pd.DataFrame([[emin, emax, min_step, max_step, mean_step, std_step]], columns=column_names)
                
                df = pd.concat([df, df2], axis = 0, ignore_index = True)
    
        if print_summary:    
            print(f'Energy Range and Energy-Step Summary for {spectra_name}')
            print(f"\tVariation in starting energy points between spectra [eV]: {df.E_min.min():.2f}-{df.E_min.max():.2f}")
            print(f"\tVariation in ending energy points between spectra [eV]: {df.E_max.min():.2f}-{df.E_max.max():.2f}")
            print(f"\tVariation in step size of energy points between spectra [eV]: {df.Mean_E_Step.min():.2f}-{df.Mean_E_Step.max():.2f}\n")
            
            if has_e0:
                print(f'Normalization Parameters for {spectra_name}')
                print(f'\tEdge Energy Range [E0]: {min(df.E0):.2f}-{max(df.E0):.2f} eV')
                print(f'\tPre-edge start/stop Suggestion [pre1/pre2]: {max([df.E_min.max()-e0,-150]):.0f}/{-50:.0f}')
                print(f'\tPost-edge start/stop Suggestion [norm1/norm2]: {75:.0f}/{min([df.E_max.min()-e0,700]):.0f}')
                print(f'\tNormalizaion order Suggestion [nnorm]: {2:.0f}')
                print(f'\tFlatten Spectra Suggestion [make_norm]: {True}\n')
            
        return df
    
    ###############################################
    
    ##### Single-Spectra Calcualtion Functions #####
    
    ###############################################
    
    def calculate_spectrum(self, scan, sample_spectra = True, ref_spectra = True):
        '''
        
        Calculate the absorption spectrum for the sample (and reference)

        Parameters
        ----------
        scan : STR
            Name of a scan to calculate the absorption spectrum of
        sample_spectra : BOOL, optional
            Flag to determine if the absorption spectrum of the sample is calcualted.
            The default is True.
        ref_spectra : BOOL, optional
            Flag to determine if the absorption spectrum of the reference is calcualted.
            The default is True.

        Returns
        -------
        None.

        '''
        
        energy_col = self.spectra[scan]['XAS Data Structure']['energy column']
        
        self.spectra[scan]['Absorption Spectra'] = {}
        
        if sample_spectra:
            # Define signals to use to calc mu sample
            sample_numerator = self.spectra[scan]['XAS Data Structure']['sample numerator']
            sample_denominator = self.spectra[scan]['XAS Data Structure']['sample denominator']
            sample_ln = self.spectra[scan]['XAS Data Structure']['sample ln']
            sample_invert = self.spectra[scan]['XAS Data Structure']['sample invert']
            
            # Extract data from signal columns
            photon_energy = self.spectra[scan]['BL Data'][energy_col]
            
            # Start update 4/17/2024
            #samp_numerator = self.spectra[scan]['BL Data'].__dict__[sample_numerator]
            #samp_denominator = self.spectra[scan]['BL Data'].__dict__[sample_denominator]
            
            # Update 4/17/2024 to include multi-channle numerators and denomunators. Check to see if there are multiple channels to be added (e.g multi-pixel Ge detetor) and sum them, otherwise use single dataset
            if isinstance(sample_numerator, list):
                for i in range(len(sample_numerator)):
                    if i == 0:
                        samp_numerator = self.spectra[scan]['BL Data'].__dict__[sample_numerator[i]]
                    else:
                        samp_numerator = np.sum([samp_numerator, self.spectra[scan]['BL Data'].__dict__[sample_numerator[i]]], axis = 0)
            else:
                samp_numerator = self.spectra[scan]['BL Data'].__dict__[sample_numerator]
            
            # Update 4/17/2024 to include multi-channle numerators and denomunators. Check to see if there are multiple channels to be added (e.g multi-pixel Ge detetor) and sum them, otherwise use single dataset
            if isinstance(sample_denominator, list):
                for i in range(len(sample_denominator)):
                    if i == 0:
                        samp_denominator = self.spectra[scan]['BL Data'].__dict__[sample_denominator[i]]
                    else:
                        samp_denominator = np.sum([samp_denominator, self.spectra[scan]['BL Data'].__dict__[sample_denominator[i]]], axis = 0)
            else:
                samp_denominator = self.spectra[scan]['BL Data'].__dict__[sample_denominator]

            # End update 4/17/2024
            
            samp_log=sample_ln
            samp_flip = sample_invert
            
            # Calcualte Absorption Coefficient 
            self.spectra[scan]['Absorption Spectra']['mu Sample'] = xfcts.create_larch_spectrum(photon_energy, samp_numerator, samp_denominator, log=samp_log, flip = samp_flip, name = scan)
        
        if ref_spectra:
            # Define signals to use to calc mu ref
            reference_numerator = self.spectra[scan]['XAS Data Structure']['reference numerator']
            reference_denominator = self.spectra[scan]['XAS Data Structure']['reference denominator']
            reference_ln = self.spectra[scan]['XAS Data Structure']['reference ln']
            reference_invert = self.spectra[scan]['XAS Data Structure']['reference invert'] 
  
            # Extract data from signal columns
            photon_energy = self.spectra[scan]['BL Data'][energy_col]
            
            # Start update 4/17/2024
            #ref_numerator = self.spectra[scan]['BL Data'].__dict__[reference_numerator]
            #ref_denominator = self.spectra[scan]['BL Data'].__dict__[reference_denominator]
            
            # Update 4/17/2024 to include multi-channle numerators and denomunators. Check to see if there are multiple channels to be added (e.g multi-pixel Ge detetor) and sum them, otherwise use single dataset
            if isinstance(reference_numerator, list):
                for i in range(len(reference_numerator)):
                    if i == 0:
                        ref_numerator = self.spectra[scan]['BL Data'].__dict__[reference_numerator[i]]
                    else:
                        ref_numerator = ref_numerator + self.spectra[scan]['BL Data'].__dict__[reference_numerator[i]]
            else:
                ref_numerator = self.spectra[scan]['BL Data'].__dict__[reference_numerator]
            
            # Update 4/17/2024 to include multi-channle numerators and denomunators. Check to see if there are multiple channels to be added (e.g multi-pixel Ge detetor) and sum them, otherwise use single dataset
            if isinstance(reference_denominator, list):
                for i in range(len(reference_denominator)):
                    if i == 0:
                        ref_denominator = self.spectra[scan]['BL Data'].__dict__[reference_denominator[i]]
                    else:
                        ref_denominator = ref_denominator + self.spectra[scan]['BL Data'].__dict__[reference_denominator[i]]
            else:
                ref_denominator = self.spectra[scan]['BL Data'].__dict__[reference_denominator]

            # End update 4/17/2024            

            ref_log = reference_ln
            ref_flip = reference_invert
            
            # Calcualte Absorption Coefficient
            self.spectra[scan]['Absorption Spectra']['mu Reference'] = xfcts.create_larch_spectrum(photon_energy, ref_numerator, ref_denominator, log=ref_log, flip = ref_flip, name = scan+'_ref')
                
        return
    
    def organize_RawData_spectrum(self, scan, remove_duplicates = True, remove_nan_inf = True, remove_zeros = True, feedback = True):
        '''
        Cleans up a spectrum's data due to high speed collection issues.
        duplicates: removes duplicate energy points caused by encoder jitter
        nan_inf: removes NAN and INF cells due to exporting/recording issues
        zeros: removes any zeros in the datafile

        Parameters
        ----------
        scan : str
            Scan Name.
        remove_duplicates : BOOl, optional
            removes duplicate values from dataset. The default is True.
        remove_nan_inf : BOOL, optional
            removes nan and inf values from dataset. The default is True.
        remove_zeros : BOOL, optional
            removes zeros from the dataset. The default is True.
        feedback : BOOL, optional
            if true, will print results to monitor. The default is True.

        Returns
        -------
        list
            Starting/ending data points prior to duplicate removal.
        list
            Starting/ending data points prior to NaN/inf removal.
        list
            Starting/ending data points prior to zero removal.

        '''
        
        energy_col = self.spectra[scan]['XAS Data Structure']['energy column']
        
        # Placeholder values if remove funciton is not called
        E_Pts1 = len(self.spectra[scan]['BL Data'][energy_col])
        E_Pts2 = E_Pts1
        naninf_Pts1 = E_Pts1
        naninf_Pts2 = E_Pts1
        Z_Pts1 = E_Pts1
        Z_Pts2 = E_Pts1
        
        
        if remove_duplicates:

            # Lenght of dateset - uses lenght of energy column
            E_Pts1 = len(self.spectra[scan]['BL Data'][energy_col])
            
            # Finds all unique values and indicies
            values, index = np.unique(self.spectra[scan]['BL Data'][energy_col], return_index=True)
            
            # Keep only index values of unique data points
            for line in self.spectra[scan]['BL Data'].array_labels:
                self.spectra[scan]['BL Data'].__dict__[line] = self.spectra[scan]['BL Data'].__dict__[line][index]
            
            # Lenght of new dataset wiht duplicates removed
            E_Pts2 = len(self.spectra[scan]['BL Data'][energy_col])
            
            if feedback:
                print(f'Duplicate data points removed for {scan}')
                print(f'Number of datapoints removed: {E_Pts1-E_Pts2}')
            

        if remove_nan_inf:
            
            print('Remove Zero still in development')
            
            
        if remove_zeros:
            
          print('Remove Zero still in development')
                
                
        return [E_Pts1, E_Pts2], [naninf_Pts1, naninf_Pts2], [Z_Pts1, Z_Pts2]
    
    
    def load_params_spectrum(self, scan, spectra_name, param_dict):
        '''
        loads paramters from dictionary using key-value into a given scans specific spectra (sample or ref) in experiment object

        Parameters
        ----------
        scan : TYPE
            DESCRIPTION.
        spectra_name : STR
            Spectra name, either 'mu Sample' or 'mu Reference'.
        param_dict : DICT
            Dictionary of key-values for paramters to be used in teh larch group for performing XAS calculations.

        Returns
        -------
        None.

        '''
        
        # Added 1/29/2024 by ASH - allows for working on a given scan in dataset, needed for real time analysis
                
        for key in param_dict.keys():
            self.spectra[scan]['Absorption Spectra'][spectra_name].__dict__[key] = param_dict[key]
            
        return
    
    def normalize_spectrum(self, scan, spectra_name):        
        '''
        Normalizes a single scan (sample or reference spectrum) in the experimentclass

        Parameters
        ----------
        scan : STR
            Name of scan to perform normalization on.
        spectra_name : STR
            "mu Sample" or "mu Reference" - which spectra in the scan should be normalized.

        Returns
        -------
        None.

        '''
        # Added 1/29/2024 by ASH - allows for working on a single scan, needed for real-time analysis
                
        xfcts.normalize_spectrum(self.spectra[scan]['Absorption Spectra'][spectra_name])

        return
    
    
    ###############################################
    
    ##### Multi-Spectra Calcualtion Functions #####
    
    ###############################################
    
    def calculate_spectra(self, sample_spectra = True, ref_spectra = True):
        '''
        Calculates the absorption spectra of every scan (and their reference channel)

        Parameters
        ----------
        scan : STR
            Name of a scan to calculate the absorption spectrum of
        sample_spectra : BOOL, optional
            Flag to determine if the absorption spectrum of the sample is calcualted.
            The default is True.
        ref_spectra : BOOL, optional
            Flag to determine if the absorption spectrum of the reference is calcualted.
            The default is True.

        Returns
        -------
        None.

        '''
        # 1/28/2024 ASH - created to try to integrate real time data analysis
        
        for key in self.spectra.keys():
            
            self.calculate_spectrum(key, sample_spectra = sample_spectra, ref_spectra = ref_spectra)
            
        return
       

    def organize_RawData(self, remove_duplicates = True, remove_nan_inf = True, remove_zeros = True, feedback = False, summary = True):
        '''
        

        Parameters
        ----------
        remove_duplicates : TYPE, optional
            DESCRIPTION. The default is True.
        remove_nan_inf : TYPE, optional
            DESCRIPTION. The default is True.
        remove_zeros : TYPE, optional
            DESCRIPTION. The default is True.
        feedback : TYPE, optional
            DESCRIPTION. The default is True.
        summary : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        '''
        
        # 1/28/2024 ASH - created to try to integrate real time data analysis 
        
        # Empty Lists for Summarry
        E_list = [] # Energy
        N_list = [] # NaN/inf
        Z_list = [] # Zeros
        
        for key in self.spectra.keys():
            
            E, N, Z = self.organize_RawData_spectrum(key, remove_duplicates = remove_duplicates, remove_nan_inf = remove_nan_inf, remove_zeros = remove_zeros, feedback = feedback)
            
            E_list.append(E)
            N_list.append(N)
            Z_list.append(Z)
        
        
        if summary:
            E_list = np.asarray(E_list)
            N_list = np.asarray(N_list)
            Z_list = np.asarray(Z_list)
            print('Data Organized')
            print(f'Range of data points per raw spectra: {E_list[:,0].min()}-{E_list[:,0].max()}')
            print(f'Range of data points after duplicates removed: {E_list[:,1].min()}-{E_list[:,1].max()}')
            print(f'Range of data points after NaN/inf removed: {N_list[:,1].min()}-{N_list[:,1].max()}')
            print(f'Range of data points after zeros removed: {Z_list[:,1].min()}-{Z_list[:,1].max()}')

        
        return
    
    
    def normalize_spectra(self, spectra_name):
        '''
        Normalzies all spectra of a type (sample/reference) in the experimental object.
        Uses normalization paramters in each scan/spectra group.

        Parameters
        ----------
        spectra_name : STR
            "mu Sample" or "mu Reference" - which spectra in the scan should be normalized.

        Returns
        -------
        None.

        '''
        # Added 1/29/2024 by ASH - allows for the single spectrum workflow fucniton
        
        for key in self.spectra.keys():
            self.normalize_spectrum(key, spectra_name)
            
        return
    
    
    def interpolate_spectra(self, start, stop, step, x_axis = 'Energy', sample = 'mu Sample'):
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
        
        # Creat list of values to interpoalte on (stop value inclusive)
        interp_E = np.arange(start, stop+step, step)
        
        # Create a df to concatinate all spectra onto
        results_df = pd.DataFrame(index = interp_E)
        results_df.index.rename(x_axis, inplace = True)
        
        for key in self.spectra.keys():
            #Write Select Data into dataframe
            time_step = self.spectra[key]['Time'] 
            
            if x_axis == 'energy':
                #y = 'flat'
                data = {x_axis:self.spectra[key]['Absorption Spectra'][sample].__dict__['energy']+self.spectra[key]['Absorption Spectra'][sample].delE,
                time_step:self.spectra[key]['Absorption Spectra'][sample].__dict__['flat']}
            
            elif x_axis == 'k':
                #y = 'chi'    
                data = {x_axis:self.spectra[key]['Absorption Spectra'][sample].__dict__['k']+self.spectra[key]['Absorption Spectra'][sample].delE,
                time_step:self.spectra[key]['Absorption Spectra'][sample].__dict__['chi']}
             
            
            temp_df = pd.DataFrame(data)
            temp_df = temp_df.set_index(x_axis)
            
            temp_interp_df = fcts.interp_df(temp_df, results_df.index)
            
            results_df = pd.concat([results_df, temp_interp_df], axis=1, join="inner")
            
        
        # Write to summary branch of DF
        self.summary[f'Interpolated {x_axis}'] = results_df
        
        return

    
    def calibrate_reference_spectra(self, edge_energy, energy_range=[-20,20], plot_range = 20, use_mean = True, overlay = True, data_filtering = True, plot_filtering = True, window_length = 5, polyorder = 2):
        '''
        Calibrates the reference channel spectra based upon edge_energy provided
        sets reference energy e0 to edge_energy
        shifts sample + reference energy scales based upon average edge energy position
        '''
        if type(energy_range) == float or type(energy_range) == int:
            energy_range = [-1*energy_range, energy_range]
        
        # Reference Lines for defining bounds of where the edge was looked for
        emin1 = edge_energy+energy_range[0]
        emax1 = edge_energy+energy_range[1]
        
        
        emin = edge_energy-plot_range
        emax = edge_energy+plot_range   
        
        e0_list = []
        delE_list = []
        groups = []
        
        # Calculate the E0 value and delE for each reference spectra
        for key in self.spectra.keys():
            # Reset DelE
            self.spectra[key]['Absorption Spectra']['mu Reference'].delE = 0.0
            self.spectra[key]['Absorption Spectra']['mu Sample'].delE = 0.0
            
            self.spectra[key]['Absorption Spectra']['mu Reference'].e0 = edge_energy
            
            # Calculate the delE from range provided in function
            # Added potential for filtering of data to find the edge for noisy/oversampled data
            e0_calc = xfcts.calculate_spectrum_e0(self.spectra[key]['Absorption Spectra']['mu Reference'], 
                                                  edge_energy, energy_range = energy_range, set_E0 = False, 
                                                  filtering = data_filtering, window_length = window_length, polyorder = polyorder)
            # Calcualte delE
            delE = edge_energy - e0_calc
            
            # Update Summary lists
            e0_list.append(e0_calc)
            delE_list.append(delE)
                
            # Store group name for plotting data
            groups.append(self.spectra[key]['Absorption Spectra']['mu Reference'])
            
            # Store Del_E if not averaging             
            if not use_mean:
                 
                self.spectra[key]['Absorption Spectra']['mu Reference'].delE = delE
                self.spectra[key]['Absorption Spectra']['mu Sample'].delE = delE
             
        # Calcualte mean delE even if not used
        del_E = edge_energy - np.asarray(e0_list).mean()
        
        # Report Energy statistics
        print('Reference Edge Finding Statistics:')
        print(f'\tReference E0 min:  {min(e0_list):.2f} eV')
        print(f'\tReference E0 max:  {max(e0_list):.2f} eV')
        print(f'\tReference E0 mean: {np.asarray(e0_list).mean():.2f} +/- {np.asarray(e0_list).std():.2f} eV')
        print(f'\tdelE range: {np.asarray(delE_list).min():.2f}-{np.asarray(delE_list).max():.2f} eV')
        print('\nMean E0 Parameters:')
        print(f'\tdelE mean: {del_E:.2f} eV')
        
        # If using mean del_E store average del_E
        if use_mean:
            for key in self.spectra.keys():  
                self.spectra[key]['Absorption Spectra']['mu Reference'].delE = delE
                self.spectra[key]['Absorption Spectra']['mu Sample'].delE = delE
        
        # Plot Data
        pfcts.plot_XANES(groups, emin, emax, spectra = 'mu', 
                         deriv = True, e0 = None, e0_line = True, ref_lines = [emin1, emax1],
                         overlay = overlay, use_legend = False, 
                         filtering = plot_filtering, window_length = window_length, polyorder = polyorder)
        
        return
        
    def find_sample_e0(self, edge_energy, energy_range = 20, use_mean = True, overlay = True, data_filtering = True, plot_filtering = True, window_length = 5, polyorder = 2):
        '''
        finds the edge positions of the samples around edge_energy value and energy_range
        set sample e0 value to found value or mean value
        use_mean = True --> use mean value of e0 for every spectra
        use_mean = False --> use calculated value of e0 for each spectra

        '''

        e0_list = [] #List to populate edge energies into
        groups = []
        
        for key in self.spectra.keys():
            e0_list.append(xfcts.calculate_spectrum_e0(self.spectra[key]['Absorption Spectra']['mu Sample'],
                                                       edge_energy, energy_range = energy_range, set_E0 = True,
                                                       filtering = data_filtering, window_length = window_length, polyorder = polyorder))
            groups.append(self.spectra[key]['Absorption Spectra']['mu Sample'])
            
        emin = edge_energy-energy_range
        emax = edge_energy+energy_range
            
        if not use_mean:
            pfcts.plot_XANES(groups, emin, emax, spectra = 'mu', 
                         deriv = True, e0 = None, e0_line = False, 
                         overlay = overlay, use_legend = False,
                         filtering = plot_filtering, window_length = window_length, polyorder = polyorder)
        
        if use_mean:
            for key in self.spectra.keys():
                self.spectra[key]['Absorption Spectra']['mu Sample'].e0 = np.asarray(e0_list).mean()
        
            print('Sample Calibraiton Statistics:')
            print(f'\tSample E0 min: {min(e0_list):.2f} eV')
            print(f'\tSample E0 max: {max(e0_list):.2f} eV')
            print(f'\tSample E0 mean: {np.asarray(e0_list).mean():.2f} +/- {np.asarray(e0_list).std():.2f} eV')
                
            pfcts.plot_XANES(groups, emin, emax, spectra = 'mu', 
                             deriv = True, e0 = None, e0_line = True, 
                             overlay = overlay, use_legend = False,
                             filtering = plot_filtering, window_length = window_length, polyorder = polyorder)
            
        return

        

    
    
    def normalize_spectra_original(self, spectra_name):        
        '''
        FILL ME IN
        spectra

        Returns
        -------
        None.

        '''
        # 1/29/2024 ASH - original function, attempting to replace with new one for real-time work flow
                
        for key in self.spectra.keys():
            xfcts.normalize_spectrum(self.spectra[key]['Absorption Spectra'][spectra_name])

        return
    
    def extract_EXAFS_spectra(self, spectra_name):
        '''
        Calculates the EXAFS spectrum from a spectrum in a larch group.
        Spectra name is either 'mu Sample' or 'mu Reference' following convention for experiment object.

        Parameters
        ----------
        spectra_name : STR
            Spectra name, either 'mu Sample' or 'mu Reference'.

        Returns
        -------
        None.

        '''        
                        
        for key in self.spectra.keys():
            xfcts.calc_spectrum_exafs(self.spectra[key]['Absorption Spectra'][spectra_name])

        return
    
    def FT_EXAFS_spectra(self, spectra_name): 
        '''
        Calculates the Fourier-Transform of the Chi data extracted from a spectrum in a larch group.
        Spectra name is either 'mu Sample' or 'mu Reference' following convention for experiment object.

        Parameters
        ----------
        spectra_name : STR
            Spectra name, either 'mu Sample' or 'mu Reference'.

        Returns
        -------
        None.

        '''
                
        for key in self.spectra.keys():
            xfcts.calc_spectrum_FT(self.spectra[key]['Absorption Spectra'][spectra_name])

        return
    
    
    def load_params(self, spectra_name, param_dict):
        '''
        loads paramters from dictionary using key-value into each spectra in experiment object

        Parameters
        ----------
        spectra_name : STR
            Spectra name, either 'mu Sample' or 'mu Reference'.
        param_dict : DICT
            Dictionary of key-values for paramters to be used in teh larch group for performing XAS calculations.

        Returns
        -------
        None.

        '''
        # 1/29/2024 ASH - updated to allow for implementaion of the scan based fucntion for real time analysis
        
        for key in self.spectra.keys():
            self.load_params_spectrum(key, spectra_name, param_dict)
            
        return

    def remove_BLData(self):
        '''
        Removes the raw data that was imported

        Returns
        -------
        None.

        '''
        for key in self.__dict__['spectra'].keys():
            self.__dict__['spectra'][key].pop('BL Data')
    
        print('Beamline data removed from experiment')
    
    #################################   
    
    ##### LCF related functions #####
    
    #################################
    
    def load_lcf_basis(self, basis_list, fit_name):
        '''
        update later basis_must be a list of larch groups with nornalized spectra
        '''
        self.analysis['LCF'][fit_name] = {}
        
        # Loads in the basis spectra
        self.analysis['LCF'][fit_name]['basis spectra'] = basis_list   
        
        return    
    
    def fit_LCF(self, fit_name, emin, emax, weights=None, minvals=0, maxvals=1, arrayname='norm', sum_to_one=False):
        '''
        TBD

        Parameters
        ----------
        xmin : TYPE
            DESCRIPTION.
        xmax : TYPE
            DESCRIPTION.
        weights : TYPE, optional
            DESCRIPTION. The default is None.
        minvals : TYPE, optional
            DESCRIPTION. The default is 0.
        maxvals : TYPE, optional
            DESCRIPTION. The default is 1.
        arrayname : TYPE, optional
            DESCRIPTION. The default is 'norm'.
        sum_to_one : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        '''
        
        minvals = [minvals]*len(self.analysis['LCF'][fit_name]['basis spectra'])
        maxvals = [maxvals]*len(self.analysis['LCF'][fit_name]['basis spectra'])
        
        self.analysis['LCF'][fit_name]['Results'] = {}
        
        for key in self.spectra.keys():
            self.analysis['LCF'][fit_name]['Results'][key] = larch.math.lincombo_fit(self.spectra[key]['Absorption Spectra']['mu Sample'],
                            self.analysis['LCF'][fit_name]['basis spectra'], weights=weights, 
                            minvals=minvals, maxvals=maxvals, arrayname=arrayname, 
                            xmin=emin, xmax=emax, sum_to_one=sum_to_one)
        
        return
   
    def lcf_report(self, fit_name):
    
        no_basis = len(self.analysis['LCF'][fit_name]['basis spectra'])

        fit_param = {'Name': [],
                    'Chi2': [],
                    'RedChi2': [],
                    'Variables': []}

        for x in range(no_basis):
            fit_param[f'Amp{x+1}'] = []
            fit_param[f'Amp{x+1}-stdev'] = []
        
        fit_param['Sum Amp'] = []
                
        for key in self.analysis['LCF'][fit_name]['Results'].keys():
            fit_param['Name'].append(key)
            fit_param['Chi2'].append(self.analysis['LCF'][fit_name]['Results'][key].chisqr)
            fit_param['RedChi2'].append(self.analysis['LCF'][fit_name]['Results'][key].redchi)
            fit_param['Variables'].append(self.analysis['LCF'][fit_name]['Results'][key].result.__dict__['nvarys'])

            # Block removed for lower suggestion given a potential bug with the appearnce of "E0_shift" in teh fit table            
            #for x, key2 in zip(list(range(no_basis)), self.analysis['LCF'][fit_name]['Results'][key].result.params.keys()):       
            #    fit_param[f'Amp{x+1}'].append(self.analysis['LCF'][fit_name]['Results'][key].result.params[key2].value)
            #    fit_param[f'Amp{x+1}-stdev'].append(self.analysis['LCF'][fit_name]['Results'][key].result.params[key2].stderr)
            
            # Jordan's Additions to account for a potential "E0_shift" term that sometimes shows up
            for x in np.arange(no_basis):
                fit_param[f'Amp{x+1}'].append(self.analysis['LCF'][fit_name]['Results'][key].result.params[f'c{x}'].value)
                fit_param[f'Amp{x+1}-stdev'].append(self.analysis['LCF'][fit_name]['Results'][key].result.params[f'c{x}'].stderr)

            
            fit_param['Sum Amp'].append(self.analysis['LCF'][fit_name]['Results'][key].result.params['total'].value)
            
                
        LCF_df = pd.DataFrame(fit_param)

        self.analysis['LCF'][fit_name]['Fit Summary'] = LCF_df
    
    
        return


    ####################################

    ##### Data Exporting Functions #####
    
    ####################################
    
    def save_normalize_spectra(self, output_directory, sep = ',', na_rep='', header = True, index = True):
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
        
        normalized_df.to_csv(output_path, sep=sep, na_rep=na_rep, header=header, index=index)
        
        return
             
    
    def save_interpXAS(self, fname_interpXAS, ext = '.csv', sep = ',', na_rep='', header = True, index = True):
        '''
        Saves the xas-process correlation dataframe to csv

        Parameters
        ----------
        fname_processparams : str
            Full path and file for saved file.

        Returns
        -------
        None.

        '''
            
        # Check file name extension
        filename, file_extension = os.path.splitext(fname_interpXAS)
        
        # Add extension if it is missing one.
        if file_extension == '':
            
                fname_interpXAS = fname_interpXAS + ext
                
        # Save the data
        self.summary['Interpolated energy'].to_csv(fname_interpXAS, sep=sep, na_rep=na_rep, header=header, index=index)
    
        print('Process Parameter Data Saved')
    
        return
    
    
    def save_processparams(self, fname_processparams):
        '''
        Saves the xas-process correlation dataframe to csv

        Parameters
        ----------
        fname_processparams : str
            Full path and file for saved file.

        Returns
        -------
        None.

        '''
            
        # Check file name extension
        filename, file_extension = os.path.splitext(fname_processparams)
        
        # Add extension if it is missing one.
        if file_extension == '':
            
                fname_processparams = fname_processparams + '.csv'
                
        # Save the data
        self.summary['XAS Spectra Process Params'].to_csv(fname_processparams, sep=',', na_rep='', header=True, index=True)
    
        print('Process Parameter Data Saved')
    
        return
    
 
    def save_lcf_results(self, fname_lcfreport, fit_name, save_spectra = True):
        '''
        Saves the LCF fit report for a given fit. If save_spectra = True, saves specta, fit, and components for each spectra in the experiment

        Parameters
        ----------
        fname_lcfreport : TYPE
            DESCRIPTION.
        fit_name : TYPE
            DESCRIPTION.
        save_spectra : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        '''
        
        # Check file name extension
        filename, file_extension = os.path.splitext(fname_lcfreport)
        
        # Add extension if it is missing one.
        if file_extension == '':
            
                fname_lcfreport = fname_lcfreport + '.csv'         
    
        # Save the LCF summary data
        self.analysis['LCF'][fit_name]['Fit Summary'].to_csv(fname_lcfreport, sep=',', na_rep='', header=True, index=True)
    
        # Saving indivisual spectra from the fit...
        if save_spectra:    
            
            # Create Subdirectory to store individual fit results               
            result_dir = os.path.dirname(fname_lcfreport)
            
            spectra_dir = fcts.create_subdir(result_dir, f'Spectra LCF Results from {fit_name}')
                
            
            for key in self.spectra.keys():
                
                # Build exportable results dataframe                
                result_dict = {
                    'Enregy': self.analysis['LCF'][fit_name]['Results'][key].xdata,
                    'Sample': self.analysis['LCF'][fit_name]['Results'][key].ydata,
                    'Fit': self.analysis['LCF'][fit_name]['Results'][key].yfit,
                    }
                
                for key2 in self.analysis['LCF'][fit_name]['Results'][key].ycomps.keys():
                    result_dict[f'component: {key2}'] = self.analysis['LCF'][fit_name]['Results'][key].ycomps[key2]
                
                df = pd.DataFrame.from_dict(result_dict)
                
                
                # Define filename of results:
                fname_spectrumLCF = key + '_LCFsummary.csv'
                
                # Define where the data will be saved
                output_path2 = os.path.join(spectra_dir, fname_spectrumLCF)
                
                # Save the LCF summary data
                df.to_csv(output_path2, sep=',', na_rep='', header=True, index=True)
                
        print('LCF Data Saved')
    
        return
        
    
    
    ### Multi - Plotting Fucntions
    
    
    def plot_XANES_spectra(self, emin, emax, samp_ref = 'mu Reference', spectra = 'mu', deriv = True, e0 = None, e0_line = True, ref_lines = None, overlay = True, use_legend = True, cmap_name = 'brg', filtering = True, window_length = 5, polyorder = 2):
        
        larch_groups = []
        
        for key in self.spectra.keys():
            larch_groups.append(self.spectra[key]['Absorption Spectra'][samp_ref])

        pfcts.plot_XANES(larch_groups, emin, emax, spectra = spectra, deriv = deriv, e0 = e0, e0_line = e0_line, ref_lines = ref_lines, overlay = overlay, use_legend = use_legend, cmap_name = cmap_name, filtering = filtering, window_length = window_length, polyorder = polyorder)
    
    
    def plot_SampRef_XANES(self, emin, emax, spectra = 'mu'):
        ### depreciated
        
        Samp_group = []
        Ref_group = []
        
        for key in self.spectra.keys():
            Samp_group.append(self.spectra[key]['Absorption Spectra']['mu Sample'])
            Ref_group.append(self.spectra[key]['Absorption Spectra']['mu Reference'])

        
        pfcts.plot_SampRef_XANES(Samp_group, Ref_group, emin, emax, spectra = spectra)
    
    
    def plot_XAS_spectra(self, emin, emax, calib_line = True): ### Plotting Data - needs improvement!!!!
        # Define Figure [2 panel side by side]
        fig1 = plt.figure(constrained_layout=True, figsize = (12,5))
        spec1 = gridspec.GridSpec(ncols = 1, nrows = 1, figure = fig1)

        f1_ax1 = fig1.add_subplot(spec1[0])
        f1_ax2 = fig1.add_subplot(spec1[0])

        #Plot Reference and Sample Spectra
        for key in self.spectra.keys():
            f1_ax1.plot(self.spectra[key]['Absorption Spectra']['mu Reference'].energy, self.spectra[key]['Absorption Spectra']['mu Reference'].mu)
            f1_ax2.plot(self.spectra[key]['Absorption Spectra']['mu Sample'].energy, self.spectra[key]['Absorption Spectra']['mu Sample'].mu)
        
            if calib_line:
                f1_ax1.plot([self.spectra[key]['Absorption Spectra']['mu Reference'].e0, 
                         self.spectra[key]['Absorption Spectra']['mu Reference'].e0],
                        [min(self.spectra[key]['Absorption Spectra']['mu Reference'].mu), max(self.spectra[key]['Absorption Spectra']['mu Reference'].mu)],
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

        
    def plot_LCF_results(self, fit_name, error_bars = True, process_parameter = None):
        # Define Figure [2 panel side by side]
        fig1 = plt.figure(constrained_layout=True, figsize = (10,10))
        spec1 = gridspec.GridSpec(ncols = 6, nrows = 10, figure = fig1)

        f1_ax1 = fig1.add_subplot(spec1[0:7,0:])
        f1_ax2 = fig1.add_subplot(spec1[7:,0:], sharex=f1_ax1)

        for col in self.analysis['LCF'][fit_name]['Fit Summary'].columns:
            if not 'stdev' in col:
                if 'Amp' in col and not 'Sum' in col:
                    if error_bars:
                        f1_ax1.errorbar(self.analysis['LCF'][fit_name]['Fit Summary'][col].index, 
                                    self.analysis['LCF'][fit_name]['Fit Summary'][col], 
                                    yerr = self.analysis['LCF'][fit_name]['Fit Summary'][col+'-stdev'],
                                    label = col)
                    else:
                        f1_ax1.plot(self.analysis['LCF'][fit_name]['Fit Summary'][col].index, 
                                    self.analysis['LCF'][fit_name]['Fit Summary'][col], 
                                    label = col)
                elif 'Sum' in col:
                    f1_ax1.plot(self.analysis['LCF'][fit_name]['Fit Summary'][col], label = col)

        f1_ax2.plot(self.analysis['LCF'][fit_name]['Fit Summary']['RedChi2'], color = 'k')


        if process_parameter != None:
            f1_ax3 = f1_ax1.twinx()
            f1_ax3.plot(self.summary['XAS Spectra Process Params'][process_parameter].values, color = 'k')

        fig1.subplots_adjust(wspace=0, hspace=0)
        
        f1_ax1.legend(loc='center right')
        f1_ax1.set_ylabel('Fraction of Basis')
        f1_ax1.set_ylim([-0.1, 1.1])
        f1_ax2.set_xlabel('Spectra')
        f1_ax2.set_ylabel('Reduced Chi^2')
        
        
    def plot_LCF(self, spectra, fit_name, emin = None, emax = None, ymin = None, ymax = None):
    
        # Plot definition
        fig1 = plt.figure(constrained_layout=True, figsize = (10,6))
        spec1 = gridspec.GridSpec(ncols = 1, nrows = 1, figure = fig1)
    
        f1_ax1 = fig1.add_subplot(spec1[0])
        
        x = self.analysis['LCF'][fit_name]['Results'][spectra].xdata
        y1 = self.analysis['LCF'][fit_name]['Results'][spectra].ydata
        y2 = self.analysis['LCF'][fit_name]['Results'][spectra].yfit
        y4 = y1-y2
    
        f1_ax1.plot(x, y1, color = 'k', label = 'Data')
        f1_ax1.plot(x, y2, color = 'r', label = 'Fit')
        f1_ax1.plot(x, y4, color = 'b', label = 'Resid')
        
    
        for key2 in self.analysis['LCF'][fit_name]['Results'][spectra].ycomps.keys():
            y3 = self.analysis['LCF'][fit_name]['Results'][spectra].ycomps[key2]
            f1_ax1.plot(x,y3, label = key2)
        
        plt.xlim(emin,emax)
        plt.ylim(ymin, ymax)
        plt.title(spectra)
        plt.legend(loc='center right', bbox_to_anchor=(1.75, 0.5))