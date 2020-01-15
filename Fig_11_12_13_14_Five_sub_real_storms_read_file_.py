"""
Created on Sun May 21 20:42:17 2017

@author: Fengwei Hung
"""
from __future__ import division
import os
os.chdir(r'C:\School\2nd Paper\model\Figure_i_adjusted_2019')  # for my computer only
import numpy as np
import matplotlib.pyplot as plt
#from cat_runoff_w_infil import sub_runoff
#from cat_runoff_w_infil import Muskingun
import pandas as pd
#from def_storm_CDO import storm_id
#from rainfall_analysis_hour import storms
from scipy.special import lambertw

#from Hydrograph_plot import Muskingun_plot 
cutoff = 0.10 # cfs
i_metric = '90perc_i'
city = 'Seattle'
#city = 'Philadelphia'
filename = city + '_Storm_dataframe_i_adjusted.xlsx'  # unit: =in/hr

tn = 3
dt =1/6
all_area = 1*43560  # surface area (ft^2) / 1 ac = 43560 ft^2
n_sub=5      # number of subwatershed
sub_p= np.array([0.2,0.2,0.2,0.2,0.2])
sub_area = all_area*sub_p
imperv_perc = 1.0   # percentage of impervious area
gi_s_p= 0.0333333           # gi surface area to subwatershed area (percentage)
gi_s = sub_area[0]*gi_s_p     # Surface area of GI (ft^2)   (1%)
DR = 6
gi_d = gi_s*DR     # drainage area (ft^2)  
Sf=[0.03,0.03,0.03]   # slope % [s_p, s_imp, s_gi]
manning_n=[0.01, 0.01, 0.01]
bern_h = 6 #in
n_yr=2013-1980+1
###########################
##  Indepent experiment
###########################
### Chaning parameters 
location = 0
L = 100            # length of the watershed

P_crit = 0.001 # in/hr
P_lag = 2
metric = 'in' # indicate the metric in the input file: mm or in
#%%  Read excel 
data_s = pd.read_excel(filename,sheetname='Data_s')
CSO_it=[]
CSO_V =[]
for i in range(n_sub+1):    
    s_name1 = 'CSO_it_%s' %i
    cso_it = pd.read_excel(filename,sheetname=s_name1)
    cso_it1 = cso_it.values
    cso_it2 = cso_it1.reshape(len(cso_it1))
    CSO_it.append(cso_it2)
    s_name2 = 'CSO_V%s' %i

    cso_v = pd.read_excel(filename,sheetname=s_name2)
    cso_v1 = cso_v.values
    cso_v2 = cso_v1.reshape(len(cso_v1))
    CSO_V.append(cso_v2)

#%% Calculate CSO reduction and the indices of the CSOs eliminated by GI
s = np.array([500,350,250,160,80,20])
location = np.arange(1,6,1)


CSO_re = []  ## index of the CSOs elimiated
CSO_V_re=[]  # the CSO volume reduction
# convert CSO_re (set) to lists
CSO_v_all = np.zeros((n_sub+1,len(CSO_V[0])))  # cso volume (no GI, GI at sub1, ...) including the CSOs eliminated (0)
CSO_v_all[0] = CSO_V[0]                        # CSO volume without GI
temp = list(CSO_it[0])                         # the index of CSO without GI
cso_i_list = []                                     # renumbering the CSO storm; CSO_it are lists of the storm id)
cso_i_list.append(range(len(CSO_it[0])))             # the index of CSOs     

for i in location:
    CSO_re.append(set(CSO_it[0])-set(CSO_it[i]))
    index_list = []                                       # the index of CSO with GI at location i
    for j in CSO_it[i]:
        index = temp.index(j)
        index_list.append(index)
    cso_i_list.append(index_list)
    CSO_v_all[i][index_list] = CSO_V[i]
    
for i in range(n_sub):
    CSO_re[i] = map(int, list(CSO_re[i]))         # CSO stroms that were fully treated
    CSO_V_re.append(CSO_v_all[0]-CSO_v_all[i+1])  # CSO reduction when GI is located at sub i, i in 1 to 5
   
#%%  Calculate the yearly ave CSO vol and the standard error
ymd = data_s.iloc[CSO_it[0]].index
#data_s.loc[data_s.index.year == 1980]
years = np.arange(1980,2014,1)
cso_y_l = [] #CSO volume for each year with GI at each location
cso_y_n = [] #CSO occurance for each year with GI at each location
# calculate the number of storms of each years

for location in range(6):
    cso_y=np.zeros(2013-1980+1)   # create a 0-array for CSO volume
    cso_n=np.zeros(2013-1980+1)   # create a 0-array for CSO occurance
    cso_ay=[]     # the CSO storms of each year given a GI location
    n0=0
    for i in range(2013-1980+1):   # loop hthough years
        storm_iy=[]                # index of storms for a particular year i
        for j in range(len(cso_i_list[location])):
            iy = CSO_it[location][j].astype(int)  # for CSO stroms j
            if years[i] == data_s.index.year[iy]:  #the year of the storm j = year i
                storm_iy.append(cso_i_list[location][j]) # stored cso id in storm_iy
        cso_ay.append(storm_iy)        #identify the storm id in scenarios with GI located at "location"
        
# calculate annual CSO volume        
        ns = len(cso_ay[i])       # the number of CSO storms of an year
        cso_y[i] = np.array(CSO_v_all[location][storm_iy]).sum()  # kgal
        cso_n[i] = ns
        n0=n0+ns
    cso_y_l.append(cso_y)
    cso_y_n.append(cso_n)
CSO_ym = np.average(cso_y_l, axis = 1)  # average annual CSO volume
CSO_se = np.std(cso_y_l,axis =1)/(n_yr**0.5) # standard error of annual CSO volume
CSO_yn = np.average(cso_y_n, axis = 1)  # no. of CSO events
CSO_n_se = np.std(cso_y_n,axis =1)/(n_yr**0.5) # standard error of annual CSO volume

h0 = []  # test if siting at S3 can provide larger SW reduction
for i in range(len(cso_y_l[0])):
    a1 = cso_y_l[3][i]<cso_y_l[5][i]
    h0.append(a1)


d = np.array([cso_y_l[0].tolist(),cso_y_l[1].tolist(),cso_y_l[2].tolist(),cso_y_l[3].tolist(),cso_y_l[4].tolist(),cso_y_l[5].tolist()])
np.savetxt(city+"_cso_y_l.csv",d, fmt='%5s',delimiter=",")
dn = np.array([cso_y_n[0].tolist(),cso_y_n[1].tolist(),cso_y_n[2].tolist(),cso_y_n[3].tolist(),cso_y_n[4].tolist(),cso_y_n[5].tolist()])
np.savetxt(city+"cso_y_n.csv",dn, fmt='%5s',delimiter=",")
# the CSO volume are highly correlated between senarios (GI siting at different locations)

#%% Figure 11 - Histogram of the CSO reduction and frequency 

CSO_yrm=[]
for i in np.arange(1,6):
    cso_y1 =cso_y_l[0]-cso_y_l[i]
    CSO_yrm.append(cso_y1*3.7854/0.40468)
    
CSO_reV = np.average(CSO_yrm, axis = 1) 
# average annual CSO volume m3/ha
CSO_V_se = 1.96*np.std(CSO_yrm  ,axis =1)/(n_yr**0.5) # standard error of annual CSO volume

CSO_yrn = []
for i in np.arange(1,6):
    cso_y1 =cso_y_n[0]-cso_y_n[i]
    CSO_yrn.append(cso_y1)
CSO_reN = np.average(CSO_yrn, axis = 1)  # average annual CSO volume    
CSO_N_se = 1.96*np.std(CSO_yrn  ,axis =1)/(n_yr**0.5) # standard error of annual CSO volume

    
    
plt.rcParams.update({'font.size': 14})
width = 0.5  
fig1= plt.figure(1, figsize=(6, 3))
plt.bar(range(5), CSO_reV, width, color='b', yerr=CSO_V_se,alpha =0.5)
plt.xlim(-0.5,5.)
plt.ylim(0,300)
#plt.yticks(np.arange(0, 300, step=50))
plt.ylabel("Volume Reduction (m^3/ha)")
plt.xlabel("GI Location")
plt.xticks( np.arange(0.,4.2,1), ('S1', 'S2', 'S3', 'S4', 'S5'), rotation='vertical' )

plt.grid()
plt.tight_layout()
plt.savefig("Fig11_%s_CSO_vol_re_dr_%s_bh_%s.png" %(city,DR,bern_h))

width = 0.5  
fig2= plt.figure(2, figsize=(6, 3))
#plt.plot(range(7), np.repeat(CSO_yn[0],7))
plt.bar(range(5), CSO_reN, width, color='b', yerr=CSO_N_se,alpha =0.5)
plt.xlim(-0.5,5)
plt.ylim(0,5)
#plt.yticks(np.arange(0, 10, step=1))

plt.ylabel("Occurrance Reduction")
plt.xlabel("GI Location")
#plt.yticks( np.arange(min(CSO_f)-1, max(CSO_f)+1, 1))
plt.xticks(  np.arange(0.,4.2,1), ('S1', 'S2', 'S3', 'S4', 'S5'), rotation='vertical' )
#plt.xticks(x, labels, rotation='vertical')
plt.grid()
plt.tight_layout()
plt.savefig("Fig11_%s_CSO_occurance_re_dr_%s_bh_%s.png" %(city,DR,bern_h))
plt.close()

#    
#%% Figure 12: GI storage utilization for CSO volume reduction
color1 =['b','y','r','c','k','m','g','r','lightgreen']
markers1 = ['o','>','<','*','x','.','1','.','|','_',' ']
N = len(data_s['v'])
XX=range(N)

for i in range(N):
    XX[i]=[data_s['length'][i], data_s[i_metric][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] 
    #  duration, characteristic i, volume, max_i, storm id
XX = np.array(XX)

plt.rcParams.update({'font.size': 14})

Alpha = 0.8
S_gi = gi_s_p*bern_h/12*7.48*43560/1000 # 5.43 kgal GI's storage

location=np.arange(1,6,1)
for loc in np.arange(0,5,1):
    cso_m = cso_i_list[loc+1]                 # index of CSOs with GI located at (loc+1)
    CSO_V_RE = CSO_V_re[loc]           # CSO reduction
    cso_a1 = np.nonzero(CSO_V_RE/S_gi<=0.4)   # index of CSOs with GI storage utilization rate <0.4
    cso_a2 = np.nonzero((CSO_V_RE/S_gi>0.4) * (CSO_V_RE/S_gi<=0.8)) # index of CSOs with GI storage utilization rate between 0.4 and 0.8
    cso_a3 = np.nonzero(CSO_V_RE/S_gi>0.8)    # index of CSOs with GI storage utilization rate >0.8
    list_a = [CSO_it[0][cso_a1], CSO_it[0][cso_a2],CSO_it[0][cso_a3]]  # lists of storm id that indicate the CSO storms in each storage utilization group
    list_b = [cso_a1, cso_a2, cso_a3]                                  # lists of CSO id that indicate the CSO reduction with GI located at (loc+1)
    label1 = ['Utilization < 40%', '40% <= Utilization <80%', '80% <= Utilization']
    plt.figure(figsize=(8, 4))  # Tr : Tn
    ax = plt.gca()
    ax.scatter(XX[CSO_it[0],2]*DR/(bern_h), XX[CSO_it[0],1]*tn*DR/(bern_h), c='white',linewidths=1, edgecolor="k",
                   marker=markers1[0], s=CSO_V[0]*10, alpha=0.4 , label = 'CSO Volume') 

    for i in np.arange(0,3,1):   
        ax.scatter(XX[list_a[i],2]*DR/(bern_h), XX[list_a[i],1]*tn*DR/(bern_h), c=color1[i],linewidths=0,
                       marker=markers1[0], s = CSO_V_RE[list_b[i]]*10, alpha=(0.4+i*0.2) , label = label1[i]) 
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.ylabel("Tn")
    plt.xlabel("Tr")
    plt.xlim(0.1,10)
    plt.ylim(0.1,10)
    plt.tight_layout()
    plt.grid()
    plt.savefig("Fig14%s_%s_TrTn_S_adjusted_i_%s.png" %((location[loc],city,location[loc])))
    plt.close()


#%% Figure 13 - use mean intensity - no Adjustment 
# CSO, nonCSO, Loc relevant, loc irrelevant storms and the CSO separation lines
color1 =['b','y','r','c','k','m','g','r','lightgreen']

Alpha = 1
# XX[ length, mean_i , volume, max_i]
N = len(data_s['v'])
XX=range(N)

for i in range(N):
    XX[i]=[data_s['length'][i], data_s['mean_i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] 
    #  duration, characteristic i, volume, max_i, storm id
XX = np.array(XX)

loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO = list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

cso_m = CSO_it[0].astype(int)

s1 = 40

qcso= cutoff
GN = tn*qcso*DR/bern_h
Gn=GN

env_y = []
Tr = np.arange(Gn+0.001,6,0.01)

for Tr1 in Tr:
    zz = -(np.exp(-Tr1/Gn)*Tr1/Gn)
#    print(zz)
    w = lambertw(zz,0)
    Tn1 = (Gn*Tr1)/(Tr1+Gn*w.real)
    env_y.append(Tn1)
#plt.plot(Tr, env_y)
# show non CSO storm events
fig3= plt.figure(3, figsize=(8,4))  # Volume vs Maximum Intensity
ax = plt.gca()
plt.plot([GN,10],GN*np.ones(2), c='b',label = 'CSO Separation Line')
plt.plot(GN*np.ones(2),[GN,10], c ='b')
plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='white',
                   marker=markers1[0], alpha=0.8, edgecolor="grey", label = 'Non-CSOs')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[1], alpha=0.9, edgecolor="k", label = 'Persisting CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',
            linewidths=0, marker=markers1[2], s=s1, alpha=1 ,edgecolor="b", label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',
            linewidths=0, marker=markers1[3], s=80, alpha =Alpha,edgecolor="r", label = 'Loc-relevant CSOs')

plt.plot(Tr,env_y,label = 'CSO Envelope',c='g')
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("Tn")
plt.xlabel("Tr")
plt.xlim(0.01,10)
plt.ylim(0.01,10)
plt.rcParams.update({'font.size': 16})
plt.legend()
plt.legend(loc="upper left", scatterpoints=1, fontsize=14)
plt.tight_layout()
plt.grid()
plt.savefig("Fig12_GN_%s_tn_%s_mean_i.png" %(city,tn))
plt.close()


#%% Figure 14-1: i*tr is replace by volume in Tr calculation
#   adjustment factor = 1, characteristic intensity = 90% intensity
Alpha = 1
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['90perc_i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# XX[ length, i , volume, max_i]
loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

cso_m = CSO_it[0].astype(int)

s1 = 40  # the circle size on the figure
qcso= cutoff
GN = tn*qcso*DR/bern_h
a = 1 # adjustion fraction for the storage increase from infiltration of long duration storms
# only appied to the horozontal bound 
env_y = []
Gn_x=GN
Gn_y=a*GN
Gn=GN
Tr = np.arange(Gn+0.001,6,0.01)

for Tr1 in Tr:
    zz = -(np.exp(-Tr1/Gn)*Tr1/Gn)
#    print(zz)
    w = lambertw(zz,0)
    Tn1 = a*(Gn*Tr1)/(Tr1+Gn*w.real)
    env_y.append(Tn1)
#plt.plot(env_y)
# show non CSO storm events
plt.rcParams.update({'font.size': 16})
fig4 = plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
ax = plt.gca()
plt.plot([Gn_x,10],(Gn_y)*np.ones(2), c='b',label = 'CSO Separation Line')   # y = HN
plt.plot(Gn_x*np.ones(2),[Gn_y,10], c ='b')                    # x = HN
plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='white',
                   marker=markers1[0], alpha=0.8, edgecolor="grey", label = 'Non-CSOs')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[1], alpha=0.9, edgecolor="k", label = 'Persisting CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[2], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',linewidths=0, marker=markers1[3], s=80, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant CSOs')
plt.plot(Tr,env_y,label = 'CSO Envelope',c='g')
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("Tn")
plt.xlabel("Tr")
plt.xlim(0.1,10)
plt.ylim(0.1,10)
#plt.legend()
#plt.legend(loc="upper left", scatterpoints=1, fontsize=14)

plt.tight_layout()

plt.grid()
plt.savefig("Fig13-1_GN_%s_tn_%s_centroid_tn_90perc.png" %(city,tn))
plt.close()
#
print'total cso no. = ', len(CSO_it[0])
print 'loc-irr = ', len(loc_irr)
print 'loc-rel = ', len(loc_rel)
print 'CSOs = ', len(CSOs)
print 'non-CSOs', len(nonCSO)
#%% Figure 14-2: i*tr is replace by volume in Tr calculation
#   adjustment factor = 0.8, characteristic intensity = 90% intensity
Alpha = 1
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['90perc_i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# XX[ length, i , volume, max_i]
loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

cso_m = CSO_it[0].astype(int)

s1 = 40  # the circle size on the figure
qcso= cutoff
GN = tn*qcso*DR/bern_h
a = 0.6 # adjustion fraction for the storage increase from infiltration of long duration storms
# only appied to the horozontal bound 
env_y = []
Gn_x=GN
Gn_y=a*GN
Gn=GN
Tr = np.arange(Gn+0.001,6,0.01)

for Tr1 in Tr:
    zz = -(np.exp(-Tr1/Gn)*Tr1/Gn)
#    print(zz)
    w = lambertw(zz,0)
    Tn1 = a*(Gn*Tr1)/(Tr1+Gn*w.real)
    env_y.append(Tn1)
#plt.plot(env_y)
# show non CSO storm events
plt.rcParams.update({'font.size': 16})
fig5 = plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
ax = plt.gca()
plt.plot([Gn_x,10],(Gn_y)*np.ones(2), c='b',label = 'CSO Separation Line')   # y = HN
plt.plot(Gn_x*np.ones(2),[Gn_y,10], c ='b')                    # x = HN
plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='white',
                   marker=markers1[0], alpha=0.8, edgecolor="grey", label = 'Non-CSOs')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[1], alpha=0.9, edgecolor="k", label = 'Persisting CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[2], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',linewidths=0, marker=markers1[3], s=80, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant CSOs')
plt.plot(Tr,env_y,label = 'CSO Envelope',c='g')
ax.set_yscale('log')
ax.set_xscale('log')

ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("Tn")
plt.xlabel("Tr")
plt.xlim(0.1,10)
plt.ylim(0.1,10)

#plt.legend()
#plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()

plt.grid()
plt.savefig("Fig13-2_GN_%s_tn_%s_Tn_Tr_dr_%s_adjusted_envelope_line_90perc.png" %(city,tn,DR))
plt.close()

#
print'total cso no. = ', len(CSO_it[0])
print 'loc-irr = ', len(loc_irr)
print 'loc-rel = ', len(loc_rel)
print 'CSOs = ', len(CSOs)
print 'non-CSOs', len(nonCSO)

