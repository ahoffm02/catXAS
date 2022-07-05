# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:09:30 2022

@author: ashoff
"""

##############################################################################

                                # Modules #
                        
##############################################################################

# File Handling
import os

#Other Functions
import re
import errno


# Data organization
import pandas as pd
import numpy as np


##############################################################################

                        # NON-XAS FUNCTIONS #
                        
##############################################################################

def create_subdir(parent_dir, sub_dir):
    '''
    Fuction to make a subdirectory in specified directory. Checks to see if 
    the directory already exists, and eitehr makes it or doesn't. Updates
    command line to the status of the subdirectory.
    
    Parameters
    ----------
    parent_dir : STR
        directory to place subdirector
    sub_dir : STR
        name of subdirectory

    Returns
    -------
    newdir : str
        full path strin of created director

    '''

    newdir = os.path.join(parent_dir, sub_dir)

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


def mergeindex(df1, df2, method = 'time'):
    """
    interpolates df2 onto df1 index
    
    df1: dataframe wiht index used for reindexing. requires datetime index
    
    df2: dataframe to be reindexed. requires datetime index
    
    Note: if df1 index exceeds limits of df2, NaN will be returned in values that can not be interpoalted
    
    Return: df2 wiht df1 index 
    """
    
    
    df2 = df2.reindex(df2.index.union(df1.index)).interpolate(method=method, limit_area = None).reindex(df1.index)
    
    return df2
