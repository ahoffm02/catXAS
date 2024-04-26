# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 16:40:38 2024

@author: ashoff
"""
# Storage for depreciated Functions until they are confimed useless


#### From Experiment Class ####

def calculate_spectra(self, sample_spectra = True, ref_spectra = True):
        '''
        docstring TBD
        '''
        #1/28/2024 ASH - Original Function - removed to sumplify the funciton using inport spectrum and another funciton that loops the funciton
        
        for key in self.spectra.keys():
            
            self.spectra[key]['Absorption Spectra'] = {}
            
            if sample_spectra:
                # Define signals to use to calc mu sample
                sample_numerator = self.spectra[key]['XAS Data Structure']['sample numerator']
                sample_denominator = self.spectra[key]['XAS Data Structure']['sample denominator']
                sample_ln = self.spectra[key]['XAS Data Structure']['sample ln']
                sample_invert = self.spectra[key]['XAS Data Structure']['sample invert']
                
                # Extract data from signal columns
                photon_energy = self.spectra[key]['BL Data'].Energy
                samp_numerator = self.spectra[key]['BL Data'].__dict__[sample_numerator]
                samp_denominator = self.spectra[key]['BL Data'].__dict__[sample_denominator]
                samp_log=sample_ln
                samp_flip = sample_invert
                
                # Calcualte Absorption Coefficient 
                self.spectra[key]['Absorption Spectra']['mu Sample'] = xfcts.create_larch_spectrum(photon_energy, samp_numerator, samp_denominator, log=samp_log, flip = samp_flip, name = key)
            
            if ref_spectra:
                # Define signals to use to calc mu ref
                reference_numerator = self.spectra[key]['XAS Data Structure']['reference numerator']
                reference_denominator = self.spectra[key]['XAS Data Structure']['reference denominator']
                reference_ln = self.spectra[key]['XAS Data Structure']['reference ln']
                reference_invert = self.spectra[key]['XAS Data Structure']['reference invert'] 
  
                # Extract data from signal columns
                photon_energy = self.spectra[key]['BL Data'].Energy
                ref_numerator = self.spectra[key]['BL Data'].__dict__[reference_numerator]
                ref_denominator = self.spectra[key]['BL Data'].__dict__[reference_denominator]
                ref_log = reference_ln
                ref_flip = reference_invert
                
                # Calcualte Absorption Coefficient
                self.spectra[key]['Absorption Spectra']['mu Reference'] = xfcts.create_larch_spectrum(photon_energy, ref_numerator, ref_denominator, log=ref_log, flip = ref_flip, name = key+'_ref')
                
        return 
    
    
def import_biologic(self, file_name):
    # Removed due to Biologic files being hard to read. Better to pre-build a dataframe and pass it to the experiment
    '''
    Loads Biologic potentilstate data into the process_params dictionary
    of the experimental class

    Parameters
    ----------
    file_name : STR of LIST
         Single string or list of full path string(s) to Biologic file(s).

    Returns
    -------
    None.

    '''

    if type(file_name) != list:
        file_name = [file_name]
    self.process_params['Biologic Data'] = process.ReadBiologicData(file_name)
    self.summary['Biologic Filename'] = file_name
    
    return

def import_spectra_data(self, xas_data_directory, xas_data_structure, print_name = False):
    '''
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
    time_stamp = xas_data_structure['time stamp']
    time_line = xas_data_structure['time on line']
    time_format = xas_data_structure['time format']
    padded = xas_data_structure['padded scan numbers']
    is_QEXAFS = xas_data_structure['is QEXAFS']
    
    
    self.summary['XAS Spectra Files'] = xfcts.CXAS_Sorted(xas_data_directory, time_stamp = time_stamp, time_line = time_line, time_format = time_format, padded = padded, is_QEXAFS = is_QEXAFS)  
    
    
    
    for index, row in self.summary['XAS Spectra Files'].iterrows():
        
        filename = row['File Name']
        file_path = row['Path']
        
        if print_name:
            print(filename)

        self.spectra[filename] = {}
    
        self.spectra[filename]['XAS Data Structure'] = xas_data_structure
        
        self.spectra[filename]['Time'] = index
        
        col_names = self.spectra[filename]['XAS Data Structure']['column names']
        
        energy_name = xas_data_structure['energy column']
        
        self.spectra[filename]['BL Data'] = read_ascii(file_path, labels=col_names)
        
        self.spectra[filename]['BL Data'].__name__ = filename
                    
        # Determine the enregy values at the start/end of the dataset
        Einit = self.spectra[filename]['BL Data'].__dict__[energy_name][0]
        Efin = self.spectra[filename]['BL Data'].__dict__[energy_name][-1]
        
        # Invert the dataset if the energy values are not increasing
        if Einit >= Efin:
            for line in col_names:
                self.spectra[filename]['BL Data'].__dict__[line] = np.flipud(self.spectra[filename]['BL Data'].__dict__[line])

        
    return


def organize_RawData_original(self, remove_duplicates = True, remove_nan_inf = True, remove_zeros = True):
    # 1/28/2024 ASH - original function, new version below
    '''
    

    Parameters
    ----------
    remove_duplicates : TYPE, optional
        DESCRIPTION. The default is True.
    remove_nan_inf : TYPE, optional
        DESCRIPTION. The default is True.
    remove_zeros : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    '''
    
    if remove_duplicates:
        E_pts = []
        for key in self.spectra.keys():
            Temp_E_Pts1 = len(self.spectra[key]['BL Data'].Energy)
            values, index = np.unique(self.spectra[key]['BL Data'].Energy, return_index=True)
            for line in self.spectra[key]['BL Data'].array_labels:
                self.spectra[key]['BL Data'].__dict__[line] = self.spectra[key]['BL Data'].__dict__[line][index]
            Temp_E_Pts2 = len(self.spectra[key]['BL Data'].Energy)
            E_pts.append([Temp_E_Pts1, Temp_E_Pts2])
                
        E_pts_arr = np.asarray(E_pts)
        print('Duplicate data points removed')
        print(f'Range of data points per raw spectra: {E_pts_arr[:,0].min()}-{E_pts_arr[:,0].max()}')
        print(f'Range of data points per duplicates removed spectra: {E_pts_arr[:,1].min()}-{E_pts_arr[:,1].max()}')


    if remove_nan_inf:
        naninf_pts = []
        for key in self.spectra.keys():
            Temp_naninf_Pts1 = len(self.spectra[key]['BL Data'].Energy)
            for line in self.spectra[key]['BL Data'].array_labels:
                f = self.spectra[key]['BL Data'].__dict__[line]
                index = [i for i, arr in enumerate(f) if not np.isfinite(arr).all()]
                for line in self.spectra[key]['BL Data'].array_labels:
                    self.spectra[key]['BL Data'].__dict__[line] = self.spectra[key]['BL Data'].__dict__[line][index]
                Temp_naninf_Pts2 = len(self.spectra[key]['BL Data'].Energy)
            naninf_pts.append([Temp_naninf_Pts1, Temp_naninf_Pts2])
        
        naninf_pts_arr = np.asarray(naninf_pts)
        print('Zero data points removed')
        print(f'Range of data points per pre-zero spectra: {naninf_pts_arr[:,0].min()}-{naninf_pts_arr[:,0].max()}')
        print(f'Range of data points per zero removed spectra: {naninf_pts_arr[:,1].min()}-{naninf_pts_arr[:,1].max()}')
    
    if remove_zeros:
        Z_pts = []
        for key in self.spectra.keys():
            Temp_Z_Pts1 = len(self.spectra[key]['BL Data'].Energy)
            for line in self.spectra[key]['BL Data'].array_labels:
                index = np.nonzero(self.spectra[key]['BL Data'].__dict__[line])
                for line in self.spectra[key]['BL Data'].array_labels:
                    self.spectra[key]['BL Data'].__dict__[line] = self.spectra[key]['BL Data'].__dict__[line][index]
                Temp_Z_Pts2 = len(self.spectra[key]['BL Data'].Energy)
            Z_pts.append([Temp_Z_Pts1, Temp_Z_Pts2])
        
        Z_pts_arr = np.asarray(Z_pts)
        print('Zero data points removed')
        print(f'Range of data points per pre-zero spectra: {Z_pts_arr[:,0].min()}-{Z_pts_arr[:,0].max()}')
        print(f'Range of data points per zero removed spectra: {Z_pts_arr[:,1].min()}-{Z_pts_arr[:,1].max()}')
        
def load_params_original(self, spectra_name, param_dict):
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
    # 1/29/2024 ASH - original function
    
    for key1 in self.spectra.keys():
        for key2 in param_dict.keys():
            self.spectra[key1]['Absorption Spectra'][spectra_name].__dict__[key2] = param_dict[key2]