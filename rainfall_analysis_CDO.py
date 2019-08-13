"""
Created on Mon Jun 27 12:15:13 2016
Defines a function that uses a precipitation time series to estimate the parameters of the rainfall generator based on Robinson and Sivapalan.  

v2_ prepare plots for supplementary materiala

@author: fengwei
"""
import matplotlib as mpl
mpl.use('agg')
import numpy as np
import matplotlib.pyplot as plt
#from histogram_storm import storm_hist
import pandas as pd

import time
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
import os
os.chdir('C:/School/2nd Paper/model')
from def_storm_CDO import storm_id


#%%Load inputs==============================================================================
start_time = time.clock()

#filename = 'Seattle_1984_20140101.csv'  # unit: kg/m2 =mm/hr
#filename = 'PhiadlephiaAirport_1980-2013.csv'  # unit: kg/m2 =mm/hr
filename = 'Philadelphia_1980_2013.csv'
city = 'Philadelphia'
#run_dates =     ['2009-01-01', '2013-12-31']
run_dates =     ['1980-01-01', '2013-12-31']
P_crit = 0.2
P_lag = 3
metric = 'in'
cutoff = 1  # minimum rainfall (mm)
df_P, data_s = storm_id(filename, run_dates, P_crit,P_lag, metric, cutoff)
result=city+'_storms_'+run_dates[0]+run_dates[1]+'.xlsx'
data_s.to_excel(result, sheet_name='sheet1', index='true')
#%%
# Histogram Storm Duration
plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8, 12))
plt.subplot(411)
data_new= data_s[data_s['v']>12.5]
#x = np.array(data_s['length'])
x = np.array(data_new['length'])

bins = np.arange(0,50,2)
x_label = 'Storm Duration (Hr)'
n,bins,patches = plt.hist(x, bins= bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
axes = plt.gca()
axes.set_xlim([0,50])
#plt.title('Storm Histogram')
plt.grid(True)
# Histogram Storm Mean Intensity
plt.subplot(412)
x = np.array(data_new['i'])
bins = np.arange(0,12,0.5)
x_label = 'Storm Mean Intensity (mm/hr)'
n,bins,patches = plt.hist(x, bins=bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
axes = plt.gca()
axes.set_xlim([0,12])
#plt.title('Storm Histogram')
plt.grid(True)
# Histogram Total Rainfall (mm)
plt.subplot(413)
x = np.array(data_new['v'])
bins=np.arange(12,70,3)
x_label = 'Total Rainfall (mm)'
n,bins,patches = plt.hist(x, bins= bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
axes = plt.gca()
axes.set_xlim([10,70])
plt.grid(True)
# Histogram Maximum Intensity (mm.hr)
plt.subplot(414)
x = np.array(data_new['max_i'])
bins = np.arange(0,25,1)
x_label = 'Maximum Intensity (mm/hr)'
n,bins,patches = plt.hist(x, bins=bins)
plt.xlabel(x_label)
plt.ylabel('Frequency')
#plt.title('Storm Histogram')
plt.grid(True)
axes = plt.gca()
axes.set_xlim([0,25])
plt.tight_layout()
plt.savefig("Fig1_Histogram_%s_%s_%s .png" %(run_dates[0],run_dates[1],filename))

#%% Cluster
#from sklearn.cluster import KMeans
#from sklearn.preprocessing import StandardScaler
##from sklearn.datasets import make_blobs
##
##X, y = make_blobs(n_samples=n_samples, random_state=random_state)
##X = data_p.values()
#N =len(data_s['i'])
#XX=range(N)
#for i in range(N):
#    XX[i]=[data_s['i'][i],data_s['v'][i],data_s['max_i'][i]] #  duration, mean i, volume, max_i
#
## Cluster
#scaler = StandardScaler()
#scaler.fit(XX)
#Y = scaler.transform(XX)
#y_pred = KMeans(n_clusters=3).fit_predict(Y)
#XX = np.array(XX)
#plt.figure(figsize=(12, 12))
#
#plt.subplot(311)
#plt.scatter(XX[:,1], XX[:,0], c=y_pred,marker='x')  # vol, mean_i
#plt.xlabel("Total Rainfall (mm)")
#plt.ylabel("Mean(i)(mm/hr)")
##plt.title("Volume vs Mean(i)")
#plt.grid()
#
#plt.subplot(312)
#plt.scatter(XX[:,2], XX[:,0], c=y_pred,marker='x')  # max_i, mean(i)
#plt.xlabel("Maximum intensity(mm/hr)")
#plt.ylabel("Mean(i)(mm/hr)")
#plt.grid()
#
#
#plt.subplot(313)
#plt.scatter(XX[:,2], XX[:,1], c=y_pred,marker='x')  # max_i, vol
#plt.ylabel("Total Rainfall (mm)")
#plt.xlabel("Maximum intensity(mm/hr)")
##plt.title("Max(i) vs Volume")
#plt.grid()


#%%

#XX=range(N)
##XX=np.array(N)
#for i in range(N):
#    XX[i]=[data_s['v'][i],data_s['max_i'][i]] #  duration, mean i, volume, max_i
#
## Cluster
#scaler = StandardScaler()
#scaler.fit(XX)
#Y = scaler.transform(XX)
#y_pred = KMeans(n_clusters=3).fit_predict(Y)
#XX = np.array(XX)
#plt.figure(figsize=(12, 12))
#
#plt.subplot(311)
#plt.scatter(XX[:,1], XX[:,0], c=y_pred,marker='x')  # vol, mean_i
#plt.xlabel("Maximum Intensity (mm/hr)")
#plt.ylabel("Total rainfall(mm)")
##plt.title("Volume vs Mean(i)")
#plt.grid()
#
#
#





