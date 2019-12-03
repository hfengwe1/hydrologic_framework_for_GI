"""
Created on Mon Jun 27 12:15:13 2016
Defines a function that uses a precipitation time series to estimate the parameters of the rainfall generator based on Robinson and Sivapalan.  
Original code provided by D. Wilusz
Plot for Figure 10
@author: Fengwei Hung
"""

import numpy as np
import pandas as pd
from matplotlib import rcParams
rcParams['font.family'] = 'serif'

#%%Load inputs==============================================================================
def storm_id(filename, run_dates, P_crit,P_lag, metric):
#    filename: precipitation data in in/hr
#    run_dates: from xxxx-xx-xx to xxxx-xx-xx
#    P_crit   #threshold to be considered part of a storm.  
#   P_lag: # the stop time between storms to be considered as separate storms 
#   metric: mm or in                
    p_crit_adj = 0.001
    df_obs = pd.read_csv(filename, index_col=0, parse_dates=True).truncate(before=run_dates[0], after=run_dates[1])    
    df_1 = df_obs[['P']].copy()#  Define input P 
    if metric == 'in':
        df_1['P'] = df_1['P']   # inch
    elif  metric == 'mm':
        df_1['P'] = df_1['P']*25.4  #mm
    else:
        print'Wrong metric!!!!'

    rng = pd.date_range(run_dates[0],run_dates[1], freq='H')  #create a time series from DATE to DATE per hour 
    ts = pd.Series(np.zeros(len(rng)), index=rng, name='P')   #create a series 'P' of 0 for the time series
    df_0 = pd.DataFrame(ts)                                   #create a data frame for the time series
    df_P =  df_0.add(df_1, fill_value=0)                      # add the values of precip. in the dataframe just created
    N_hr = len(df_P)                                          #calculate the number of the total hours
#==============================================================================
#%%Identify storm events 
#==============================================================================
#storm event criteria: minimum dry period of p_crit hours
#P_crit=2.0 #split into two storms if longer than this #hours.  
    df_P['storm_id']=np.zeros(N_hr)
    df_P['break_id']=np.zeros(N_hr)
    storm_count=0
    break_count=0
#first hour
    if (df_P['P'][0]>=P_crit):
        storm_count+=1
        df_P['storm_id'][0]=storm_count
    else: 
        break_count+=1
        df_P['break_id'][0]=break_count    
    
    for i in np.arange(1,N_hr):
        if (df_P['P'][i]>=P_crit):
            #check if new storm or not
            if (df_P['storm_id'][i-1]==0):                
                #new storm
                storm_count+=1
            #update storm count
            df_P['storm_id'][i]=storm_count
        else:
            #check if new break
            if (df_P['break_id'][i-1]==0):
                #new break
                break_count+=1
            #update break count
            df_P['break_id'][i]=break_count
#%%  Combine storms within s_break hour
    index_storm = np.unique(df_P['storm_id'], return_index=True)[1][1:]
    index_break = np.unique(df_P['break_id'], return_index=True)[1][1:]   
    break_i = df_P['break_id'][index_break]
    storm_i = df_P['storm_id'][index_storm]

    break_length = [df_P[df_P['break_id']==x].count()[0]  for x in break_i]
    storm_length = [df_P[df_P['storm_id']==x].count()[0]  for x in storm_i]

    if (df_P['P'][0]>=P_crit):      # the first hour is a storm
        for i in np.arange(0,len(break_length)):
            if break_length[i] <=P_lag:
                s_id = df_P['storm_id'][index_storm[i]]
                df_P['storm_id'][index_break[i]:index_break[i+1]] = s_id
                df_P['break_id'][index_break[i]:index_break[i+1]] = 0
    else:                           # the first hour is a break
        for i in np.arange(1,len(break_length)):
            if break_length[i] <=P_lag:
                s_id = df_P['storm_id'][index_storm[i-1]]
                df_P['storm_id'][index_break[i]:index_break[i+1]] = s_id
                df_P['break_id'][index_break[i]:index_break[i+1]] = 0

    index_storm = np.unique(df_P['storm_id'], return_index=True)[1][1:]
    index_break = np.unique(df_P['break_id'], return_index=True)[1][1:]   
    break_i = df_P['break_id'][index_break]
    break_length = [df_P[df_P['break_id']==x].count()[0]  for x in break_i]
#

#==============================================================================
#%% compile storm breaks, duration into new dataframe
#==============================================================================   
    storm_1 = df_P['storm_id'][index_storm].tolist()
    break_1 = df_P['break_id'][index_break].tolist()
#    print(list(set(storm_1)-set(break_1)))
    for i in range(len(storm_1)):   # reassign storm and break id
#        print('New Storm_ID = %s' %(i+1))
#        print('Old Storm_ID = %s' %storm_1[i])
        df_P['storm_id'][df_P['storm_id']== storm_1[i]] = i+1
        df_P['break_id'][df_P['break_id']== break_1[i]] = i+1
    if len(break_1) > len(storm_1):
        df_P['break_id'][df_P['break_id']==break_1[len(storm_1)]] = len(storm_1)+1    
    storm_2 = df_P['storm_id'][index_storm].tolist()

# storm statistics
    max_i = []
    storm_i = []
    P80i = []
    P90i= []
    storm_v = []
    storm_length = []
    
    for x in storm_2: 
        storm_des_p = df_P['P'][df_P['storm_id']==x].describe([.8,.9])
        adj_duration =  (df_P['P'][df_P['storm_id']==x]>= p_crit_adj).sum() # adjust the duration by removing p < 0.02
        if adj_duration == 0:
            storm_length.append(1)
        else:
            storm_length.append(adj_duration)
        storm_i.append(storm_des_p['mean'])
        P80i.append(storm_des_p['80%'])
        P90i.append(storm_des_p['90%'])
        storm_v.append(df_P[df_P['storm_id']==x]['P'].sum())
        max_i.append(storm_des_p['max'])
                     
    data_s={'storm_id':np.arange(len(index_storm))+1,
            'length':storm_length,
            'mean_i':storm_i,
            'v':storm_v,
            '80perc_i':P80i,
            '90perc_i':P90i,
            'max_i':max_i}
    df_storm=pd.DataFrame(data=data_s, index=df_P.index[index_storm])
    return df_P, df_storm