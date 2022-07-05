# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:17:58 2022

@author: ashoff
"""
##############################################################################

                                # Modules #
                        
##############################################################################

# File Handling
import os
import glob2

# Timestamps
import datetime
from datetime import datetime as dt

# Data organization
import pandas as pd
import numpy as np

# X-ray Science
import larch

#From Catxas
import general as fcts

##############################################################################

                        # XAS GENERIC FUNCTIONS #
                        
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

def create_larch_spectrum(photon_energy, numerator, denominator, log=True, flip = False, name = None):
    '''
    Geenrates a X-ray Larch Group populated with an energy array (energy), and an absoprtion coefficinet array (mu).
    
    Mu is calcuated based upon list inputs that are passed to the calc_mu function.

    Parameters
    ----------
    photon_energy : LIST
        List of int/float corresponding to the photon energy of the spectrum.
    numerator : LIST
        List of int/float used as the numerator term to calcualte mu. 
        Passed to calc_mu function.
    denominator : LIST
        List uof int/float used as the denominator term to calcualte mu. 
        Passed to calc_mu function.
    log : BOOL, optional
        Bool used to determine if the natural log is applies to the calculation 
        of mu. Passed to calc_mu function. 
        The default is True.
    flip : BOOL, optional
        Bool used to determine if the numerator and denominator need to be 
        flipped duringthe calculation of mu. Passed to calc_mu function. 
        The default is False.
    name : STRING, optional
        Name of Larch group - WARNING - use of the same name in a code will 
        overwrite parameters, best left as default. The default is None.

    Returns
    -------
    spectrum : Larch Group
        Larch group with energy and mu variables, Group name is either self 
        or user assinged.

    '''
               
    if name != None:
        spectrum = larch.Group(name = name)
    elif name == None:
        spectrum = larch.Group()
    
    
    spectrum.energy = photon_energy
    spectrum.mu = calc_mu(numerator, denominator, log=True, flip = False)
    spectrum.delE = 0.0
    
    Einit = spectrum.energy[0]
    Efin = spectrum.energy[-1]
            
    if Einit >= Efin:
        for key in spectrum.__dict__.keys():
            if key not in ['__name__', 'delE']:
                spectrum.__dict__[key] = np.flipud(spectrum.__dict__[key])
    
    #if name != None:
    #    spectrum.name = name
    
    return spectrum

def calculate_spectrum_e0(larch_group, edge_energy, energy_range = 20, set_E0 = True):
    '''
    Finds the "Edge position" defined by the maximum in the first derivative.
    Note: No data smoothing applied when finding the value.
    
    Sets the E0 value in the larch gorup based up on the found value.
    
    Returns the edge value.

    Parameters
    ----------
    larch_group : Larch Group
        Larch gorup containing "energy", and "mu" arrays.
    edge_energy : float/int
        Approximate estiamted edge energy. Used to reduce dataset to a more reasonable search area.
    energy_range : float/int, optional
        +/- ernergy arange around edge_enery to look for a maximum in the dervative. The default is 20.

    Returns
    -------
    e0_pos : float
        Edge enregy found and stored in the larch group.

    '''
    
    emin = edge_energy-energy_range
    emax = edge_energy+energy_range
    
    energy = larch_group.energy + larch_group.delE # Adds any energy shift to the system
    
    ind_min = fcts.find_nearest(energy, emin)[0]
    ind_max = fcts.find_nearest(energy, emax)[0]

    e0_pos = larch.xafs.find_e0(energy[ind_min:ind_max], 
                       mu = larch_group.mu[ind_min:ind_max],
                       group = None)

    if set_E0:
        larch_group.e0 = e0_pos

    return e0_pos

def CXAS_Sorted(files_directory, time_stamp = True, time_line = 0, time_format = '%m/%d/%Y %I:%M:%S %p', padded = True, is_QEXAFS = False):
    '''
    ### GENERAL STATEMENT HERE ###

    Parameters
    ----------
    files_directory : STR
        PATH TO DIRECTORY CONTAINING ONLY XAS SPECTRA FILES.
    time_stamp : BOOL, optional
        TRUE/FALSE STATEMENT IF THE XAS SPECTRA FILES CONTAIN A TIMESTAMP IN 
        THE HEADER. The default is True.
    time_line : INT, optional
        LINE IN HEADER WHICH CONTAINTS THE TIMESTAMP WHEN THE XAS SPECTRA 
        WAS COLLECTED, ONLY USED WHEN time_stamp = True. The default is 0.
    time_format : STR, optional
        FORMAT OF THE LINE IN THE HEADER CONTAINING THE TIMESTAMP, FOLLOWS 
        DATETIME FORMAT STYLE, ONLY USED WHEN time_stamp = True. 
        The default is '%m/%d/%Y %I:%M:%S %p'.
    padded : BOOL, optional
        TRUE/FALSE STATEMENT IF THE XAS SPECTRA FILENAMES IN files_directory
        CONTAIN SPECTRA NUMBERS THAT ARE PADDED WITH ZEROS TO KEEP A CONSTANT
        CHARACTER COUNT. The default is True.

    Returns
    -------
    TOS : TYPE
        DESCRIPTION.

    '''
    
    # Use glob2 to get a list of all files in files_directory
    files = glob2.glob(files_directory + '/*')
    
    
    path_series = pd.Series(files)
    
    filename_list = []
    padded_list = []
    time_list = []
    
    for line in files:
        
        filename = os.path.basename(line)[:-4]
        filename_list.append(filename)
        
        if not padded:
            scan_no = fcts.get_trailing_number(filename)
            char = len(str(scan_no))
            filename = filename[:-char] + str(scan_no).zfill(4)
            
        padded_list.append(filename)
        
        # Open the file in read only mode
        with open(line, 'r') as f:
            count = 0
            # Read all lines in the file one by one
            for line in f:
                # For each line, check if line contains the string
                if not is_QEXAFS:
                    if count == time_line:
                        #print(line)
                        date = dt.strptime(line, time_format)
                        #print(date)
                        time_list.append(date)
                        break
                    else:
                        count = count + 1
                elif is_QEXAFS:
                    if count == time_line:
                        #print(line)
                        date = dt.strptime(line, time_format)
                        #print(date)
                        #time_list.append(date)
                    if not line.startswith('#'):
                        micros = int(line.split(sep = '\t')[-2])
                        
                        micros_to_date = datetime.timedelta(microseconds = micros)
                        
                        scan_time = date + micros_to_date
                        
                        time_list.append(scan_time)
                        
                        break
                    else:
                        count = count + 1

    filename_series = pd.Series(filename_list)
    padded_series = pd.Series(padded_list)
    
    if time_stamp:
        time_series = pd.Series(time_list)
        elapsed_time = time_series - time_series[0]
        
        temp_dict = {'Time': time_series, 'TOS [s]': elapsed_time.dt.total_seconds(), 'File Name': filename_series, 'Padded Name': padded_series, 'Path': path_series}
        
        TOS = pd.DataFrame(temp_dict)
        TOS.sort_values(by = 'Time', ignore_index = True, inplace = True)
        TOS.set_index('Time', inplace = True)
    
    else:
        time_series = pd.Series(time_list)
        
        temp_dict = {
            'Time': time_series,
            'File Name': filename_series,
            'Padded Name': padded_series,
            'Path': path_series
            }
        
        TOS = pd.DataFrame(temp_dict)
        TOS.sort_values(by = 'Padded Name', ignore_index = True, inplace = True)
        TOS.index.rename('Scan', inplace = True)
        
    return TOS

def Summarize_Energy_Range(larch_group, print_summary = False):
    '''
    larch_group
    
    return
    energy_summary = dictionary of parameters definign the range and steps of the dataset
    '''
    energy = larch_group.energy + larch_group.delE #accounts for delE varaition in data
    
    energy_summary = {'E_min': 0.0, 'E_max': 0.0, 
                      'Min_E_Step': 0.0, 'Max_E_Step': 0.0, 
                      'Mean_E_Step': 0.0, 'STD_E_Step': 0.0}
    
    energy_summary['E_min'] = min(energy)
    energy_summary['E_max'] = max(energy)
    
    temp_step = []
    
    for i in range(len(energy)-1):
        temp_step.append(energy[i+1]-energy[i])
        
    energy_summary['Min_E_Step'] = min(temp_step)
    energy_summary['Max_E_Step'] = max(temp_step)
    energy_summary['Mean_E_Step'] = np.asarray(temp_step).mean()
    energy_summary['STD_E_Step'] = np.asarray(temp_step).std()
    
    if print_summary:
        print(f'Starting Energy Value [eV]:\t{energy_summary["E_min"]:.2f}')
        print(f'Ending Energy Value [eV]:\t{energy_summary["E_max"]:.2f}')
        print(f'Smallest Energy Step [eV]:\t{energy_summary["Min_E_Step"]:.2f}')
        print(f'Largest Energy Step [eV]:\t{energy_summary["Max_E_Step"]:.2f}')
        print(f'Mean Energy Step [eV]:\t{energy_summary["Mean_E_Step"]:.2f}')
        print(f'Deviation in Energy Step [eV]:\t{energy_summary["SDT_E_Step"]:.2f}')
        
        
        print('\n\t'+'\u0332'.join("Reccomended Normalization Settings:"))
        print(f'\tpre1: {max([-150, energy_summary["E_min"]-larch_group.e0])}')
        print(f'\tpre2: {-50}')
        print(f'\tnorm1: {75}')
        print(f'\tnorm2: {min([700, round((energy_summary["E_max"]-larch_group.e0)/10)*10])}')
        print(f'\tnnorm: {2}')
        print(f'\tmake_flat: {True}')
        
    return energy_summary 

def update_norm_params(larch_group, pre1 = -100, pre2 = -50, norm1 = 75, norm2 = 300, nnorm = 2, make_flat = True):
    larch_group.pre1 = pre1
    larch_group.pre2 = pre2
    larch_group.norm1 = norm1
    larch_group.norm2 = norm2
    larch_group.nnorm = nnorm
    larch_group.make_flat = make_flat
    
def normalize_spectrum(larch_group):
    try:
        energy = larch_group.energy + larch_group.delE
        mu = larch_group.mu
        group = larch_group
        e0 = larch_group.e0
        pre1 = larch_group.pre1
        pre2 = larch_group.pre2
        norm1 = larch_group.norm1
        norm2 = larch_group.norm2
        nnorm = larch_group.nnorm
        make_flat = larch_group.make_flat
        
        larch.xafs.pre_edge(energy, mu = mu, group = group, e0 = e0, 
                                pre1 = pre1, pre2 = pre2, 
                                norm1 = norm1, norm2 = norm2, 
                                nnorm = nnorm, make_flat = make_flat)
    except:
        print('One of the normalzation parameters (pre1, pre2, norm1, norm2, nnorm, or make_flat) is not defined')

def update_autobk_params(larch_group, rbkg = 1, nknots = None, kmin = 0, kmax = None, kweight = 1, dk = 0.1,
                         win = 'hanning', nfft = 2048, kstep = 0.05, k_std = None, chi_std = None, nclamp = 2,
                        clamp_lo = 1, clamp_hi = 1, err_sigma = 1):
    larch_group.rbkg = rbkg
    larch_group.nknots = nknots
    larch_group.kmin = kmin
    larch_group.kmax = kmax
    larch_group.kweight = kweight
    larch_group.dk = dk
    larch_group.win = win
    larch_group.nfft = nfft
    larch_group.kstep = kstep
    larch_group.k_std = k_std
    larch_group.chi_std = chi_std
    larch_group.nclamp = nclamp
    larch_group.clamp_lo = clamp_lo
    larch_group.clamp_hi = clamp_hi
    larch_group.err_sigma = err_sigma
    
def calc_spectrum_exafs(larch_group):
    try:
        energy = larch_group.energy + larch_group.delE
        mu = larch_group.mu 
        group = larch_group
        rbkg = larch_group.rbkg
        e0 = larch_group.e0
        
        edge_step = larch_group.edge_step
        nknots = larch_group.nknots
        kmin = larch_group.kmin
        kmax = larch_group.kmax
        kweight = larch_group.kweight
        
        dk = larch_group.dk
        win = larch_group.win
        nfft = larch_group.nfft
        kstep = larch_group.kstep
        k_std = larch_group.k_std
        
        chi_std = larch_group.chi_std
        nclamp = larch_group.nclamp
        clamp_lo = larch_group.clamp_lo
        clamp_hi = larch_group.clamp_hi
        
        err_sigma = larch_group.err_sigma
        
        larch.xafs.autobk(energy, mu, group = group, rbkg = rbkg, e0 = e0, 
                          edge_step = edge_step, nknots = nknots, kmin = kmin, kmax = kmax, kweight = kweight, 
                          dk = dk, win = win, nfft = nfft, kstep = kstep, k_std = k_std,
                          chi_std = chi_std, nclamp = nclamp, clamp_lo = clamp_lo, clamp_hi = clamp_hi,
                          err_sigma = err_sigma)
    except:
        print('One of the autobk parameters is incorrectly define [see Larch documentation]')

def update_xftf_params(larch_group, rmax_out = 10, kweight = 2, kmin = 3, kmax = 12, dk = 5, dk2 = 5, 
                           window = 'haning', nfft = 2048, kstep = 0.05):
    larch_group.rmax_out = rmax_out
    larch_group.kweight = kweight
    larch_group.kmin = kmin
    larch_group.kmax = kmax
    larch_group.dk = dk
    larch_group.dk2 = dk2
    larch_group.window = window
    larch_group.nfft = nfft
    larch_group.kstep = kstep

def calc_spectrum_FT(larch_group):
    k = larch_group.k
    chi = larch_group.chi
    group = larch_group
    rmax_out = larch_group.rmax_out
    kweight = larch_group.kweight
    kmin = larch_group.kmin
    kmax = larch_group.kmax
    dk = larch_group.dk
    dk2 = larch_group.dk2
    window = larch_group.window
    nfft = larch_group.nfft
    kstep = larch_group.kstep

    larch.xafs.xftf(k, chi = chi, group = group, rmax_out = rmax_out, kweight = kweight, kmin = kmin, kmax = kmax,
        dk = dk, dk2 = dk2, window = window, nfft = nfft, kstep = kstep)