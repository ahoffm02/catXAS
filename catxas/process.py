# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 11:20:13 2022

@author: ashoff
"""
##############################################################################

                                # Modules #
                        
##############################################################################

#Other Functions
import re
import os

# Timestamps
from datetime import datetime as dt

# Data organization
import pandas as pd

# Data Parsing
from itertools import islice


##############################################################################

                        # Process Stream FUNCTIONS #
                        
##############################################################################

def Read_Hiden_QGAPro(filename):
    """    
    Reads a XLSX file supplied from the  Hiden Quant Software.
    
    filename: string to file containing mass spectrometer data
    
    Returns: Pandas dataframe with a datetime index and signals matched with signal value (m/z) 
    """
    
    # generate filename without path and extension
    base_name = os.path.basename(filename)
    fname_no_ext = os.path.splitext(base_name)[0]
    
    # searching filename for time stamp
    match_str = re.search('\d+-\d+-\d+ \d+-\d+-\d+ \w+', fname_no_ext)
  
    # Convert timestamp string into datetime object
    time_zero = dt.strptime(match_str.group(), '%m-%d-%Y %I-%M-%S %p')#.date()
    
    
    # Get Metadata Header info and store it in list
    header = []
    count = 0
    with open(filename, 'r') as myfile:
        for line in myfile:
            count += 1
            header.append(line.strip())
            if line.strip() == 'data':
                #count += 1
                break
  
    # Close file
    myfile.close()
        
    # Create Dataframe after skipping header lines
    raw_data = pd.read_csv(filename, sep = ',', skiprows = count, header = [0,1])#, index_col = 0)

    
    # Fill in missing Headers for multiindex
    for i in range(len(list(raw_data))):
        if 'Elapsed' in list(raw_data)[i][0]:
            raw_data.rename(columns = {list(raw_data)[i][0]:''}, inplace = True, level=0, errors='raise')
            raw_data.rename(columns = {list(raw_data)[i][1]:''}, inplace = True, level=1, errors='raise')
        elif "Unnamed:" in list(raw_data)[i][0]:
            raw_data.rename(columns = {list(raw_data)[i][0]:list(raw_data)[i-1][0]}, inplace = True, level=0, errors='raise')
            
    # Set the elapsed time to index
    raw_data.set_index('', inplace = True)
    raw_data.rename_axis('Time', axis = 0, inplace = True)    
    
        
    # use time_zero and elapsed time to redefine the index as datatime opbkect
    
    new_index = pd.Series(pd.to_timedelta(list(raw_data.index), unit='s'), name = 'Time')

    for i in new_index.index:
        new_index.iloc[i] = time_zero + new_index.iloc[i]
        #new_index.iloc[i] = new_index.iloc[i].strftime('%m/%d/%Y %I:%M:%S %p')
    
    new_index.rename('Time', inplace = True)
    
    raw_data.index._data = new_index.values
    
    raw_data.index = pd.to_datetime(raw_data.index)
    
    return raw_data

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
    
    # Remove duplicate index
    data = data.loc[~data.index.duplicated(), :]
    
    return data

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
    
    #Clean up the string setpoints for the valves
    # 'A' = 1
    # 'B' = 2
    # 'Moving...' = 0
    
    text_columns = ['SV1 SP','SV1 Feedback','SV2 SP', 'SV2 Feedback']
    
    for line in text_columns:
        if line in data.columns:
            data[line].replace(to_replace = {'A':1, 'B':2,'A - Left position':1, 'B - Center position':2, 'Moving or Left':0}, inplace = True)
    
    # Remove duplicate index
    data = data.loc[~data.index.duplicated(), :]
    
    return data

def ReadBiologicData(filename):
    """
    tbd
    filename type = list
    """
    if isinstance(filename, list) == False:
        filename = [filename]
    # Extract Data into pandas dataframe
    for i in range(len(filename)):
        if i == 0:
            #print(filename[i])
            data = pd.read_csv(filename[i], sep = '\t', encoding= 'unicode_escape')
            data.drop(data.columns[data.columns.str.contains('unnamed',case = False)],
              axis = 1, inplace = True)
    
            # Convert "Date and Time" to timedelta + date and time file was created
            data.rename(columns={'time/s':'Time'}, inplace=True)
            data['Time'] = pd.to_datetime(data['Time'])
    
            # Reindex the data to the time stamp     
            data = data.set_index('Time')
        
        else:
            #print(filename[i])
            data2 = pd.read_csv(filename[i], sep = '\t', encoding= 'unicode_escape')
            data2.drop(data2.columns[data2.columns.str.contains('unnamed',case = False)],
              axis = 1, inplace = True)
    
            # Convert "Date and Time" to timedelta + date and time file was created
            data2.rename(columns={'time/s':'Time'}, inplace=True)
            data2['Time'] = pd.to_datetime(data2['Time'])
    
            # Reindex the data to the time stamp     
            data2 = data2.set_index('Time')
            
            data = pd.concat([data,data2], axis=0, join="inner")
            
        
        # Sort data by Index
        data.sort_index(axis=0, level=None, ascending=True, inplace=True, kind='quicksort', na_position='last', sort_remaining=True, ignore_index=False, key=None)
    
    return data