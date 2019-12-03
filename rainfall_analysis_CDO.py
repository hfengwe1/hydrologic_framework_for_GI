"""
Created on Mon Jun 27 12:15:13 2016
Defines a function that uses a precipitation time series to estimate the parameters of the rainfall generator based on Robinson and Sivapalan.  

v2_ prepare plots for supplementary materiala

@author: fengwei
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
import os
os.chdir('C:/School/2nd Paper/model')  # change the working directory
from def_storm_CDO import storm_id
start_time = time.clock()

#%%Load inputs==============================================================================
    # City name can only be "Philadelphia" or "Seattle"; this could be expanded when more data sets are added.
    # Data set need to be placed in the same folder
    #    filename: precipitation data in in/hr
    #    metric: output unit - mm or in       
    #    run_dates: from xxxx-xx-xx to xxxx-xx-xx
    #    P_crit   # threshold to be considered part of a storm in output metric 
    #    P_lag:  the stop time between storms to be considered as separate storms (hour)
city = 'Philadelphia'        # city name
#city = 'Seattle'

if city == 'Philadelphia':   
    filename = 'Philadelphia_1980_2013.csv' # data unit in/hr
elif city == 'Seattle':
    filename = 'Seattle_1980_2013.csv'  # unit: # unit: in/hr
else:
    print('precipitation data does not exist')

metric = 'mm'   # mm or in         
run_dates =     ['1980-01-01', '2013-12-31']   # data start date and end date    
P_crit = 0.2    # minimum precipitation for a storm
P_lag = 3       # hour        

df_P, data_s = storm_id(filename, run_dates, P_crit,P_lag, metric)  # generate a dataframe of the storm characteristics
    # df_P: precipitation record
    # data_s: a dataframe constains columns of 
    # max_i: max intensity, 80perc_i: 80% intensity, 90perc_i: 90% intensity, mean_i: mean intensity, 
    # storm_id: storm id, length: storm duration, v: total volume of the storm
result = city+'_storms_'+run_dates[0]+run_dates[1]+'.xlsx'  # save the storm characteristics dataframe to excel
data_s.to_excel(result, sheet_name='sheet1', index='true')
#data_s = pd.read_excel(open(result, 'rb'))   # retrive the results by reading the excel file
#%%
# storms with total volume less than 5 mm is removed from this analysis
data_new= data_s[data_s['v']>5]  

# Histograms
plt.rcParams.update({'font.size': 18})
plt.figure(figsize=(6, 12))
    # Histogram Storm Mean Intensity
plt.subplot(411)
x = np.array(data_new['length']) #
bins = np.arange(0,50,2)
x_label = 'Storm Duration (h)'
n,bins,patches = plt.hist(x, bins= bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
axes = plt.gca()
axes.set_ylim([0,1000])  # set y-axis
axes.set_xlim([0,50])   # set x-axis
axes.set_yticks(np.arange(0, 1000, step=250))
axes.set_xticks(np.arange(0, 50, step=10))
plt.grid(True)

    # Histogram Total Rainfall (mm)
plt.subplot(412)
x = np.array(data_new['v'])
bins=np.arange(5,70,3)
x_label = 'Total Rainfall (mm)'
n,bins,patches = plt.hist(x, bins= bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
axes = plt.gca()
axes.set_ylim([0,600])  # set y-axis
axes.set_xlim([5,70])   # set x-axis
axes.set_yticks(np.arange(0, 610, step=200))
axes.set_xticks(np.arange(10, 70, step=10))
plt.grid(True)

    # Histogram Storm Mean Intensity
plt.subplot(413)
x = np.array(data_new['mean_i'])
bins = np.arange(0,12,0.5)
x_label = 'Storm Mean Intensity (mm/h)'
n,bins,patches = plt.hist(x, bins=bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
axes = plt.gca()
axes.set_ylim([0,850])  # set y-axis
axes.set_xlim([0,12])   # set x-axis
axes.set_yticks(np.arange(0, 850, step=200))
axes.set_xticks(np.arange(0, 12, step=2))
plt.grid(True)

    # Histogram Maximum Intensity (mm.hr)
plt.subplot(414)
x = np.array(data_new['max_i'])
bins = np.arange(0,25,1)
x_label = 'Maximum Intensity (mm/h)'
n,bins,patches = plt.hist(x, bins=bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
plt.grid(True)
axes = plt.gca()
axes.set_ylim([0,750])  # set y-axis
axes.set_xlim([0,25])   # set x-axis
axes.set_yticks(np.arange(0, 780, step=250))
axes.set_xticks(np.arange(0, 25, step=5))
plt.tight_layout()
plt.savefig("Fig1_Histogram_%s_%s_%s_5mm .png" %(run_dates[0],run_dates[1],filename))







