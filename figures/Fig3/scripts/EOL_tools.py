'''
This is a set of functions used to process the Tecan Sunrise
data, to calculate area under the curve and compare AUCs.
'''

import pandas as pd
from itertools import product
import numpy as np
import datetime
import numpy as np


def read_rawdata(file, strain, plate):
    '''
    Reads output table from tecan sunrise and gets it into
    "tidy data" mode. Adds columns with strain and plate information
    '''
    # read dataframe
    df = pd.read_csv(file,
                     sep='\t',
                     skiprows=1,
                     index_col=0)
    
    # drop last column if nan 
    df = df.dropna(axis='columns')
    
    # melt into tidy data format
    df = pd.melt(df, ignore_index=False)
    
    # reset index to make wells a column
    df = df.reset_index()
    
    # convert string into seconds variable
    df['variable'] = df['variable'].str.strip('s')
    
    # rename columns approp
    df = df.rename(columns = {'variable':'seconds',
                              'index':'well',
                              'value':'OD600'})
    
    # convert to timedelta
    df['timedelta'] = pd.to_timedelta(df['seconds'].astype(int), unit='s')
    
    # clean columns
    df = df[['well', 'timedelta', 'OD600']]

    # add strain and plate data
    df['strain'] = [strain] * len(df)
    df['plate'] = [plate] * len(df)
    
    return df


def read_experiment(file, strain, plate, metadata_file, clean=True, crop_min=600):
    '''
    read table with read_rawdata, merge with metadata
    and clean up
    
    '''
    # => READ
    # read experiment data
    df_data = read_rawdata(file, strain, plate)
    
    # read original table
    df_metadata = pd.read_csv(metadata_file,
                              sep='\t')
    
    # merge
    df_exp = pd.merge(df_data,
                      df_metadata,
                      on=['well', 'plate'])
    
    # => CLEAN
    if clean:
        # remove non-control wells with no phage
        df_exp = df_exp[df_exp['phage'].notnull()]
        # keep only wells with culture, (remove phage neg ctrls)
        df_exp = df_exp[df_exp['culture'] == 1]
    
    if crop_min:
        # crop after specified time
        df_exp = df_exp[df_exp['timedelta'] < datetime.timedelta(minutes=crop_min)]
    
    df_exp = df_exp[['well', 'strain', 'phage', 'timedelta', 'OD600']]
    df_exp = df_exp.reset_index(drop=True)
    
    return df_exp
        
def auc_experiment(df_experiment, mean=True):
    '''
    Calculates area under the curve for each well,
    from an experiment dataframe.
    It can return the auc per well,
    or the mean/sem per phage-strain pair 
    '''
    auc_results = []

    for well in df_experiment['well'].unique():
        # subset per well, sort
        df_well = df_experiment[df_experiment['well'] == well]
        df_well = df_well.sort_values('minutes')

        # ground by substracting lowest value
        min_val = np.array(df_well['OD600'].to_list()).min()
        df_well['OD600_norm'] = df_well['OD600'] - min_val

        # convert to list
        od_list = df_well['OD600_norm'].to_list()

        # do the average for two consecutive entries
        sigma = []

        for i in range(len(od_list) - 1):
            partial_avg = (od_list[i + 1] + od_list[i]) / 2
            sigma.append(partial_avg)

        # make last value 0
        sigma.append(0)

        # sum all values
        auc = np.array(sigma).sum()

        # append data to results
        phage = df_well['phage'].unique()[0]
        strain = df_well['strain'].unique()[0]

        auc_results.append([well, phage, strain, auc])
    
    df_auc = pd.DataFrame(auc_results)
    df_auc.columns = ['well', 'phage', 'strain', 'auc']
    
    if mean:
        # group and get mean
        df_auc_mean = df_auc.groupby(['strain', 'phage']).agg({'auc': 'mean'}).reset_index()
        # group and get sem
        df_auc_mean['auc_sem'] = df_auc.groupby(['strain', 'phage']).agg({'auc': 'sem'}).reset_index()['auc']
        # rename columns
        df_auc_mean.columns = ['strain', 'phage', 'auc_mean', 'auc_sem']
        
        return df_auc_mean
    
    else:
        return df_auc

    
def las_experiment(df_experiment):
    '''
    Calculate liquid assay score, by substracting the
    experiment auc from the control auc, dividing by the
    control auc and mult. by 100
    '''
    # get mean area under the curve for experiment
    df_auc = auc_experiment(df_experiment)
    
    las_results = []
    
    # liquid assay score calculated per strain
    for strain in df_auc['strain'].unique():
        # subset to strain
        df_strain = df_auc[df_auc['strain'] == strain].copy()
        
        # extract control growth value
        control_auc = df_strain[df_strain['phage'] == '-']['auc_mean'].to_list()[0]
        
        # calculate liquid assay score
        df_strain['las'] =  ((control_auc - df_strain['auc_mean']) / control_auc) * 100
        las_results.append(df_strain)
    
    df_las = pd.concat(las_results)
    
    return df_las
