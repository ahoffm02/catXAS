# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 17:53:53 2022

@author: ashoff
"""

##############################################################################

                                # Modules #
                        
##############################################################################

# Data organization
import numpy as np

# Plotting
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

# Data Processing
from scipy.signal import savgol_filter

#From Catxas
import general as fcts


##############################################################################

                        # Generic Plotting FUNCTIONS #
                        
##############################################################################

def get_cmap(n, name='brg'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def plot_XYs_Vline(xy_list, vline_pos, emin, emax, axis_label = ['Photon Energy (eV)','Norm mux'], size = [12, 10]):
    
    # Colors:
    cmap = get_cmap(len(xy_list), name = 'Set1')    
    
    # Set up plot   
    fig1 = plt.figure(constrained_layout=True, figsize = (size[0], size[1]))
    spec1 = gridspec.GridSpec(ncols = 1, nrows = 1, figure = fig1)
    ax1 = fig1.add_subplot(spec1[0,0])
    
    for i in range(len(xy_list)):
        
        x = xy_list[i][0]
        y = xy_list[i][1]
        
        ax1.plot(x, y, color = cmap(i), linestyle = 'solid')
        
    ax1.axvline(vline_pos, color = 'k')
    
    # Set Plot Formatting
    ax1.set_xlim([emin, emax])
    ax1.set_xlabel(axis_label[0])
    ax1.set_ylabel(axis_label[1])
    
    return    



##############################################################################

                        # XAS Plotting FUNCTIONS #
                        
##############################################################################


def plot_XANES(larch_groups, emin, emax, spectra = 'mu', deriv = True, e0 = None, e0_line = True, ref_lines = None, overlay = True, use_legend = True, cmap_name = 'brg', filtering = True, window_length = 5, polyorder = 2):
    '''
    UPDATE!!!!
    
    Plot XANES Spectra from Larch Grouops

    Parameters
    ----------
    larch_groups : LIST
        Larch group or list of Larch groups to be plotted.
    emin : FLAOT
        Starting energy to plot from.
    emax : FLOAT
        Endign energy to plot to.
    spectra : STR, optional
        Type of spectra to plot. Common options are 'mu', 'flat', 'norm'. The default is 'mu'.
    deriv : BOOL, optional
        Option to include the derivatie plot next to main spectra. The default is True.
    e0 : FLOAT or None, optional
        E0 value to use when adding the E0 line. None will look for E0 in the larch group. If larch group does nto have e0 value, e0_line is forced to False. The default is None.
    e0_line : BOOL, optional
        Option to include a line vertically located at the E0 value specified in the function or from the Larch group e0 field. The default is True.

    Returns
    -------
    None.

    '''
    # Store value of e0 from input
    e0_temp = e0
    
    
    # Make the group a list if it isnt already one
    if type(larch_groups) != list:
        larch_groups = [larch_groups]
    
    # Count number of groups in list    
    num_groups = len(larch_groups)
    
    # Make the rel_lines a list if it isnt already one
    if ref_lines != None and type(ref_lines) != list:
        ref_lines = [ref_lines]
    
    # Set Figure paramters based on number of groups and use of derivative plot
    if deriv:
        # 2 plots wide
        width = 12
        ncols = 2
        
    elif not deriv:
        # 1 plots wide
        width = 6
        ncols = 1
        
    if overlay:
        # 1 plots tall
        height = 5
        nrows = 1
        
    
    elif not overlay:
        # n plots tall
        height = num_groups*5
        nrows = num_groups
        
    # Set Color Grid
    cmap = get_cmap(num_groups, name = cmap_name)
        
    # Set up plot   
    fig1 = plt.figure(constrained_layout=True, figsize = (width, height))
    spec1 = gridspec.GridSpec(ncols = ncols, nrows = nrows, figure = fig1)
    
    # Add subplots and dta for each group
    for i in range(num_groups):
        # Dignostic check, print line number
        #print(i)
        
        # Define X (energy) and Y (absorption coefficient) 
        x1 = larch_groups[i].energy + larch_groups[i].delE
        
        if filtering:
            y1 = savgol_filter(larch_groups[i].__dict__[spectra], window_length = window_length, polyorder = polyorder) 
        elif not filtering: 
            y1 = larch_groups[i].__dict__[spectra]
        
        # Calcualte Der. - will be smoothed already if filtering applied above
        x2 = (x1[:-1] + x1[1:]) / 2
        der_y = np.diff(y1) / np.diff(x1)
        
                
        # Check to See if E0 has been defined in spectra if e0 == None. If not, ignore plotting
        if e0_line and e0_temp == None:
            try:
                e0 = larch_groups[i].e0
                # Diagnostic check, print found e0 value from group
                #print(e0)
            except:
                e0 = None
                print(f'E0 not defined in {larch_groups[i].__name__} or function - E0 not plotted')
                # Diagnostic check, alert if there was no E0 foune in group
                #print('no E0 value found')
        
        # Plotting  
        if overlay:
        
            # Update plotting limits with each spectra
            if i == 0:
                # Build the subplot
                ax1 = fig1.add_subplot(spec1[0,0])
                
                # Determine the current window ranges for X and Y axis
                E_range = [fcts.find_nearest(x1, emin)[0], fcts.find_nearest(x1, emax)[0]]
                mu_min = min(y1[E_range[0]:E_range[1]])
                mu_max = max(y1[E_range[0]:E_range[1]])
                del_mu = abs(mu_max-mu_min)
                
            else:
                
                # Determine the current window ranges for X and Y axis
                E_range = [fcts.find_nearest(x1, emin)[0], fcts.find_nearest(x1, emax)[0]]
                temp_mu_min = min(y1[E_range[0]:E_range[1]])
                temp_mu_max = max(y1[E_range[0]:E_range[1]])
                
                # Determine if the current window ranges need to be updated
                mu_min = min(mu_min, temp_mu_min)
                mu_max = max(mu_max, temp_mu_max)
                del_mu = abs(mu_max-mu_min)
                
            # Add Spectra to Subplot
            ax1.plot(x1, y1, label = larch_groups[i].__name__, color = cmap(i), linestyle = 'solid')
            
            # Add E0 line if appropriate
            if e0_line and e0 != None:
                #print('adding e0 line')
                ax1.axvline(e0, color = 'k')
                #ax1.plot([e0, e0], [mu_min, mu_max], color = cmap(i))
                
            # Add reference lines if appropriate
            if ref_lines != None:
                for line in ref_lines:
                    ax1.axvline(line, color = 'r')
                    #ax1.plot([line, line], [mu_min, mu_max], color = 'r')
                
            if use_legend:
                ax1.legend(loc="upper right")
            
            # Set Plot Formatting
            ax1.set_xlim([emin, emax])
            ax1.set_xlabel('Photon Energy (eV)')

            ax1.set_ylim([mu_min-0.25*del_mu, mu_max+0.25*del_mu])
            ax1.set_ylabel(spectra)
            
            # Add Deriv Subplot + Data
            if deriv:
                
                # Update plotting limits with each spectra
                if i == 0:
                    
                    # Create Subplot for Derivative
                    ax2 = fig1.add_subplot(spec1[0,1])
                    
                    # Determine the window size                    
                    der_range = [fcts.find_nearest(x2, emin)[0], fcts.find_nearest(x2, emax)[0]]
                    der_mu_min = min(der_y[der_range[0]:der_range[1]])
                    der_mu_max = max(der_y[der_range[0]:der_range[1]])
                    der_del_mu = abs(der_mu_min-der_mu_max)

                else:
                    
                    # Determine the window size
                    der_range = [fcts.find_nearest(x2, emin)[0], fcts.find_nearest(x2, emax)[0]]
                    temp_der_mu_min = min(der_y[der_range[0]:der_range[1]])
                    temp_der_mu_max = max(der_y[der_range[0]:der_range[1]])
                    
                    # Determine if window size needs to be updated
                    der_mu_min = min(der_mu_min, temp_der_mu_min)
                    der_mu_max = max(der_mu_max, temp_der_mu_max)
                    der_del_mu = abs(der_mu_max-der_mu_min)
                
                # Add Data to Axis
                ax2.plot(x2, der_y, label = larch_groups[i].__name__, color = cmap(i), linestyle = 'solid')
                
                if use_legend:
                    ax2.legend(loc="upper left")
                
                # Set Plot Formatting
                ax2.set_xlim([emin, emax])
                ax2.set_xlabel('Photon Energy (eV)')

                ax2.set_ylim([der_mu_min-0.25*der_del_mu, der_mu_max+0.25*der_del_mu])  
                ax2.set_ylabel(f'd{spectra}/dE')
                
                # Add E0 line if appropriate
                if e0_line and e0 != None:
                    ax2.axvline(e0, color = 'k')
                    #ax2.plot([e0, e0],[der_mu_min-0.25*abs(der_del_mu), der_mu_max+0.25*abs(der_del_mu)],color = 'k')
                
                # Add reference lines if appropriate
                if ref_lines != None:
                    for line in ref_lines:
                        ax2.axvline(line, color = 'r')
                        #ax2.plot([line, line], [der_mu_min-0.25*abs(der_del_mu), der_mu_max+0.25*abs(der_del_mu)], color = 'r')
            
        if not overlay:
            
            # Create Subplot for Spectra
            ax1 = fig1.add_subplot(spec1[i,0])
            
            #Find Limits to Adjsut Axis Scales
            E_range = [fcts.find_nearest(x1, emin)[0], fcts.find_nearest(x1, emax)[0]]
            mu_min = min(y1[E_range[0]:E_range[1]])
            mu_max = max(y1[E_range[0]:E_range[1]])
            del_mu = abs(mu_max-mu_min)
            
            # Add Spectra to Subplot
            ax1.plot(x1, y1, label = larch_groups[i].__name__, color = cmap(i), linestyle = 'solid')
            
            if use_legend:
                ax1.legend(loc="upper right")
            
            # Set Plot Formatting
            ax1.set_xlim([emin, emax])
            ax1.set_xlabel('Photon Energy (eV)')

            ax1.set_ylim([mu_min-0.25*del_mu, mu_max+0.25*del_mu])
            ax1.set_ylabel(spectra)
            
            # Add E0 line if appropriate
            if e0_line and e0 != None:
                ax1.axvline(e0, color = 'k')
                #ax1.plot([e0, e0], [mu_min, mu_max], color = 'k')
                
            # Add reference lines if appropriate
            if ref_lines != None:
                for line in ref_lines:
                    ax1.axvline(line, color = 'r')
                    #ax1.plot([line, line], [mu_min, mu_max], color = 'r')
        
            # Add Deriv Subplot + Data
            if deriv:
                # Create Subplot for Derivative
                ax2 = fig1.add_subplot(spec1[i,1])
                
                
                #Find Limits to Adjsut Axis Scales
                der_range = [fcts.find_nearest(x2, emin)[0], fcts.find_nearest(x2, emax)[0]]
                der_mu_min = min(der_y[der_range[0]:der_range[1]])
                der_mu_max = max(der_y[der_range[0]:der_range[1]])
                der_del_mu = abs(der_mu_min-der_mu_max)
                    
                # Add Data to Axis
                ax2.plot(x2, der_y, label = larch_groups[i].__name__, color = cmap(i), linestyle = 'solid')
                
                if use_legend:
                    ax2.legend(loc="upper left")
                
                # Set Plot Formatting
                ax2.set_xlim([emin, emax])
                ax2.set_xlabel('Photon Energy (eV)')

                ax2.set_ylim([der_mu_min-0.25*der_del_mu, der_mu_max+0.25*der_del_mu])  
                ax2.set_ylabel(f'd{spectra}/dE')
                
                # Add E0 line if appropriate
                if e0_line and e0 != None:
                    ax2.axvline(e0, color = 'k')
                    #ax2.plot([e0, e0],[der_mu_min-0.25*abs(der_del_mu), der_mu_max+0.25*abs(der_del_mu)],color = 'k')
                
                # Add reference lines if appropriate
                if ref_lines != None:
                    for line in ref_lines:
                        ax2.axvline(line, color = 'r')
                        #ax2.plot([line, line], [der_mu_min-0.25*abs(der_del_mu), der_mu_max+0.25*abs(der_del_mu)], color = 'r')
                        
    return
    
def plot_chi(larch_groups, kweight = 2, kmin = 0, kmax = 15, overlay = True, use_legend = True, cmap_name = 'brg'):
    
    # Make the group a list if it isnt already one
    if type(larch_groups) != list:
        larch_groups = [larch_groups]
    
    # Count number of groups in list    
    num_groups = len(larch_groups)
    
    # Set Color Grid
    cmap = get_cmap(num_groups, name = cmap_name)
    
    # Set Figure paramters based on number of groups and other inputs
    if overlay:
        # 1 plot wide, 1 plot tall
        width = 6
        ncols = 1
        
        # n plots tall
        height = 5
        nrows = 1
        
    elif not overlay:
        # 1 plots wide
        width = 6
        ncols = 1
        
        # n plots tall
        height = num_groups*5
        nrows = num_groups
        
    # Set up plot   
    fig1 = plt.figure(constrained_layout=True, figsize = (width, height))
    spec1 = gridspec.GridSpec(ncols = ncols, nrows = nrows, figure = fig1)
    
    # Add subplots and data for each group
    for i in range(num_groups):
        
        # Define X (k) and Y (chi*k^n) 
        k = larch_groups[i].k
        chi = np.multiply(larch_groups[i].chi, k**kweight)
        
        
        # Overlaying Data:
        if overlay:
            
            # Update plotting limits with each spectra
            if i == 0:
                
                # Add Spectra Subplots + Data
                ax1 = fig1.add_subplot(spec1[0,0])
                
                # Set Limits                
                k_range = [fcts.find_nearest(k, kmin)[0], fcts.find_nearest(k, kmax)[0]]
                chi_min = min(chi[k_range[0]:k_range[1]])
                chi_max = max(chi[k_range[0]:k_range[1]])
                
            else:
                k_range = [fcts.find_nearest(k, kmin)[0], fcts.find_nearest(k, kmax)[0]]
                temp_chi_min = min(chi[k_range[0]:k_range[1]])
                temp_chi_max = max(chi[k_range[0]:k_range[1]])
            
                chi_min = min(chi_min, temp_chi_min)
                chi_max = max(chi_max, temp_chi_max)
            
            # Add Data to Axis
            ax1.plot(k, chi, label = larch_groups[i].__name__, color = cmap(i), linestyle = 'solid')
            
            if use_legend:
                ax1.legend(loc="upper left")
            
            # Add Formatting
            ax1.set_xlim([kmin, kmax])
            ax1.set_xlabel('k [A^-1]')

            ax1.set_ylim([1.25*chi_min, 1.25*chi_max])
            ax1.set_ylabel(f'k^{kweight} weighted X(k)')

        # Overlaying Data:
        if not overlay:
            
            # Add Spectra Subplots + Data
            ax1 = fig1.add_subplot(spec1[i,0])
            
            # Define X (k) and Y (chi*k^n) 
            k = larch_groups[i].k
            chi = np.multiply(larch_groups[i].chi, k**kweight)
            
            # Plotting limits with each spectra
            k_range = [fcts.find_nearest(k, kmin)[0], fcts.find_nearest(k, kmax)[0]]
            chi_min = min(chi[k_range[0]:k_range[1]])
            chi_max = max(chi[k_range[0]:k_range[1]])
            
            # Add Data to Axis
            ax1.plot(k, chi, label = larch_groups[i].__name__, color = cmap(i), linestyle = 'solid')
            
            if use_legend:
                ax1.legend(loc="upper left")
            
            # Add Formatting
            ax1.set_xlim([kmin, kmax])
            ax1.set_xlabel('k [A^-1]')

            ax1.set_ylim([1.25*chi_min, 1.25*chi_max])
            ax1.set_ylabel(f'k^{kweight} weighted X(k)')

    
def plot_FT(larch_groups, Rmin = 0, Rmax = 6, magnitude = True, imaginary = True, real = False, overlay = True, use_legend = True, cmap_name = 'brg'):
    
    # Make the group a list if it isnt already one
    if type(larch_groups) != list:
        larch_groups = [larch_groups]
    
    # Count number of groups in list    
    num_groups = len(larch_groups)
    
    # Set Color Grid
    cmap = get_cmap(num_groups, name = cmap_name)
    
    # Set Figure paramters based on number of groups and other inputs
    if overlay:
        # 1 plot wide, 1 plot tall
        width = 6
        ncols = 1
        
        # n plots tall
        height = 5
        nrows = 1
   
    elif not overlay:
        # 1 plots wide
        width = 6
        ncols = 1
        
        # n plots tall
        height = num_groups*5
        nrows = num_groups
    
    # Set up plot   
    fig1 = plt.figure(constrained_layout=True, figsize = (width, height))
    spec1 = gridspec.GridSpec(ncols = ncols, nrows = nrows, figure = fig1)
    
    # Add subplots and data for each group
    for i in range(num_groups):
        
        # Define X (R) and Y (magnitude, imaginary, and/or real) 
        R = larch_groups[i].r
        if magnitude:
            mag =  larch_groups[i].chir_mag
        if imaginary:
            imag = larch_groups[i].chir_im
        if real:
            rmag = larch_groups[i].chir_re
        
        # Overlaying Data:
        if overlay:
             
            # Set subplot and Update plotting limits with each spectra
            if i == 0:
                
                # Add Spectra Subplots + Data
                ax1 = fig1.add_subplot(spec1[0,0])
                
                R_Range = [fcts.find_nearest(R, Rmin)[0], fcts.find_nearest(R, Rmax)[0]]
                FT_max = max(mag[R_Range[0]:R_Range[1]])
                
            else:
                R_Range = [fcts.find_nearest(R, Rmin)[0], fcts.find_nearest(R, Rmax)[0]]
                temp_FT_max = max(mag[R_Range[0]:R_Range[1]])
            
                FT_max = max(FT_max, temp_FT_max)
                
             # Add Data to Axis
            if magnitude:
                ax1.plot(R, mag, label = str(larch_groups[i].__name__ + '-magnitude'), color = cmap(i), linestyle = 'solid')
            if imaginary:
                ax1.plot(R, imag, label = str(larch_groups[i].__name__ + '-imaginary'), color = cmap(i), linestyle = 'dashed')
            if real:
                ax1.plot(R, rmag, label = str(larch_groups[i].__name__ + '-real'), color = cmap(i), linestyle = 'dotted')
            
            if use_legend:
                ax1.legend(loc="upper right")
            
            
            # Add Formatting
            ax1.set_xlim([Rmin, Rmax])
            ax1.set_xlabel('R [A]')

            ax1.set_ylim([-1.25*FT_max, 1.25*FT_max])
            ax1.set_ylabel(f'k^{larch_groups[i].kweight} weighted X(R)')
        
        # Overlaying Data:
        if not overlay:
            
            # Add Spectra Subplots + Data
            ax1 = fig1.add_subplot(spec1[i,0])
            
            # Plotting limits with each spectra
            R_Range = [fcts.find_nearest(R, Rmin)[0], fcts.find_nearest(R, Rmax)[0]]
            FT_max = max(mag[R_Range[0]:R_Range[1]])
            
            # Add Data to Axis
            if magnitude:
                ax1.plot(R, mag, label = str(larch_groups[i].__name__ + '-magnitude'), c = 'k')
            if imaginary:
                ax1.plot(R, imag, label = str(larch_groups[i].__name__ + '-imaginary'), c = 'b')
            if real:
                ax1.plot(R, rmag, label = str(larch_groups[i].__name__ + '-real'), c = 'r')
            
            if use_legend:
                ax1.legend(loc="upper right")
            
            
            # Add Formatting
            ax1.set_xlim([Rmin, Rmax])
            ax1.set_xlabel('R [A]')

            ax1.set_ylim([-1.25*FT_max, 1.25*FT_max])
            ax1.set_ylabel(f'k^{larch_groups[i].kweight} weighted X(R)')
