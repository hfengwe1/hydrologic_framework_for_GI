"""
Created on Sun May 21 20:42:17 2017

@author: Fengwei Hung
"""
from __future__ import division
import os
os.chdir(r'C:\School\2nd Paper\model\Figure_SI_0730_NEW\Fig 0414 Envolope Lines\Real storm_philly')  # for my computer only
import numpy as np
import matplotlib.pyplot as plt
#from cat_runoff_w_infil import sub_runoff
#from cat_runoff_w_infil import Muskingun
import pandas as pd
#from def_storm_CDO import storm_id
#from rainfall_analysis_hour import storms
from scipy.special import lambertw

#from Hydrograph_plot import Muskingun_plot 
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
cutoff = 0.20 # cfs
bern_h = 8 #in
n_yr=2013-1980+1
###########################
##  Indepent experiment
###########################
### Chaning parameters 
location = 0
L = 100            # length of the watershed
city = 'Philladelphia'
filename = 'Philadelphia_Storm_dataframe_S_8in.xlsx'  # unit: =in/hr

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

#%%    
markers1 = ['o','>','<','1','x','.','1','*','|','_',' ']
s = np.array([500,350,250,160,80,20])
location = np.arange(1,6,1)

N = len(data_s['i'])
XX=range(N)
CSO_re = []  ## CSO event reduction
CSO_V_re=[]
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# convert CSO_re (set) to lists
for i in location:
    CSO_re.append(set(CSO_it[0])-set(CSO_it[i]))
for i in range(n_sub):
    CSO_re[i] = map(int, list(CSO_re[i]))  # CSO stroms that were fully treated
    CSO_V_re.append(CSO_V[0]-CSO_V[i+1])
    
#%% plot no 3: focus on manageable storms
color1 =['b','y','r','y','c','k','m','g','r','lightgreen']

Alpha = 0.8
S_gi = gi_s_p*bern_h/12*7.48*43560/1000 # 5.43 kgal GI's storage
# remove extreme storms

## GI at S1   
location=np.arange(1,6,1)
#cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
cso_m = CSO_it[0].astype(int)
#cso_i_m = np.arange(len(CSO_it[0]))  # the index of the storms satisfying the criteria

plt.rcParams.update({'font.size': 14})

for loc in np.arange(0,5,1):

    CSO_V_RE = CSO_V_re[loc][cso_m]
    cso_a1 = np.nonzero(CSO_V_RE/S_gi<=0.4)
    cso_a2 = np.nonzero((CSO_V_RE/S_gi>0.4) * (CSO_V_RE/S_gi<=0.8))
    cso_a3 = np.nonzero(CSO_V_RE/S_gi>0.8)
    cso_c1 = cso_m[cso_a1].astype(int)
    cso_c2 = cso_m[cso_a2].astype(int)
    cso_c3 = cso_m[cso_a3].astype(int)
    list_a = [cso_a1, cso_a2,cso_a3]
    list_c = [cso_c1,cso_c2,cso_c3]
    label1 = ['Utilization < 40%', '40% <= Utilization <80%', '80% <= Utilization']
    plt.figure(figsize=(8, 4))  # Tr : Tn
    ax = plt.gca()
    ax.scatter(XX[cso_m,2]*DR/(bern_h), XX[cso_m,1]*tn*DR/(bern_h), c='grey',linewidths=0,
                   marker=markers1[0], s=CSO_V[0][cso_m]*10, alpha=Alpha , label = 'CSO Volume') 

    for i in np.arange(0,3,1):
    
        ax.scatter(XX[list_c[i],2]*DR/(bern_h), XX[list_c[i],1]*tn*DR/(bern_h), c=color1[i],linewidths=0,
                       marker=markers1[0], s=CSO_V_RE[list_a[i]]*10, alpha=Alpha , label = label1[i]) 
    ax.set_yscale('log')
    ax.set_xscale('log')
#    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    plt.xticks(, ('No GI','S1', 'S2', 'S3', 'S4', 'S5'), rotation='vertical' )

    plt.ylabel("Tn")
    plt.xlabel("Tr")
    plt.xlim(0.1,10)
    plt.ylim(0.01,10)
#    plt.legend(loc="lower left", markerscale=1., scatterpoints=1, fontsize=10)   
#    lgnd = plt.legend(loc="upper left", scatterpoints=1, fontsize=13)
#    lgnd.legendHandles[0]._sizes = [50]
#    lgnd.legendHandles[1]._sizes = [50]
#    lgnd.legendHandles[2]._sizes = [50]
#    lgnd.legendHandles[3]._sizes = [50]   

 #change the marker size manually for both lines

    plt.tight_layout()
    plt.grid()
    plt.savefig("Fig2%s_%s_TrTn_S_%s.png" %((location[loc],city,location[loc])))



#%%  Calculate the yearly ave CSO vol and the standard error
ymd = data_s.iloc[CSO_it[0]].index
data_s.loc[data_s.index.year == 1980]

years = np.arange(1980,2014,1)
cso_y=np.zeros(2013-1980+1)
cso_y_l = [] #CSO volume for each year with GI at each location
cso_y_n = [] #CSO occurance for each year with GI at each location

# calculate the number of storms of each years

for location in range(6):
    cso_y=np.zeros(2013-1980+1)
    cso_n=np.zeros(2013-1980+1)
    storm_ay=[]     # the CSO storms of each year given a GI location
    n0=0
    for i in range(2013-1980+1):   # loop hthough years
        storm_iy=[]                # index of storms for a particular year i
        for j in CSO_it[location].astype(int):  # for CSO stroms j
            if years[i] == data_s.index.year[j]:  #the year of the storm j = year i
                storm_iy.append(j)   # stored storm j in storm_iy
        storm_ay.append(storm_iy)        #identify the storm # in scenarios that GI sitted at each subwatershed
        
# calculate annual CSO volume        
        ns = len(storm_ay[i])       # the number of CSO storms of an year
        cso_y[i] = np.array(CSO_V[location][storm_iy]).sum()  # kgal
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
np.savetxt("cso_y_l.csv",d, fmt='%5s',delimiter=",")
dn = np.array([cso_y_n[0].tolist(),cso_y_n[1].tolist(),cso_y_n[2].tolist(),cso_y_n[3].tolist(),cso_y_n[4].tolist(),cso_y_n[5].tolist()])
np.savetxt("cso_y_n.csv",dn, fmt='%5s',delimiter=",")
# the CSO volume are highly correlated between senarios (GI siting at different locations)

#%% Plot CSO reduction

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
fig1= plt.figure(4, figsize=(6, 3))
plt.bar(range(5), CSO_reV, width, color='b', yerr=CSO_V_se,alpha =0.5)
plt.xlim(-0.5,5.)
plt.ylim(0,250)
plt.yticks(np.arange(0, 300, step=50))
plt.ylabel("Volume Reduction (m^3/ha)")
plt.xlabel("GI Location")
plt.xticks( np.arange(0.,4.2,1), ('S1', 'S2', 'S3', 'S4', 'S5'), rotation='vertical' )

plt.grid()
plt.tight_layout()
plt.savefig("Fig1_%s_CSO_vol_re_dr_%s_bh_%s.png" %(city,DR,bern_h))

width = 0.5  
fig1= plt.figure(1, figsize=(6, 3))
#plt.plot(range(7), np.repeat(CSO_yn[0],7))
plt.bar(range(5), CSO_reN, width, color='b', yerr=CSO_N_se,alpha =0.5)
plt.xlim(-0.5,5)
plt.ylim(0,5)
plt.yticks(np.arange(0, 6, step=1))

plt.ylabel("Occurrance Reduction")
plt.xlabel("GI Location")
#plt.yticks( np.arange(min(CSO_f)-1, max(CSO_f)+1, 1))
plt.xticks(  np.arange(0.,4.2,1), ('S1', 'S2', 'S3', 'S4', 'S5'), rotation='vertical' )
#plt.xticks(x, labels, rotation='vertical')
plt.grid()
plt.tight_layout()
plt.savefig("Fig1_%s_CSO_occurance_re_dr_%s_bh_%s.png" %(city,DR,bern_h))



#%% length vs mean i

color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
markers1 = ['o','>','<','1','x','.','1','*','|','_',' ']
s = np.array([500,350,250,160,80,20])
N = len(data_s['i'])
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# convert CSO_re (set) to lists
for i in range(n_sub):
    CSO_re[i] = map(int, list(CSO_re[i]))

fig2= plt.figure(2, figsize=(6, 3))

plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(12, 12))
location = np.arange(1,6,1)
Alpha =  1 # transparebcy

plt.subplot(311)

plt.scatter(XX[CSO_it[0].astype(int),0], XX[CSO_it[0].astype(int),1], c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO')  # vol, mean_i
for i in range(n_sub):
    plt.scatter(XX[CSO_re[i],0], XX[CSO_re[i],1], c=color1[i],linewidths=0,
                   marker=markers1[0], s=s[i], alpha=Alpha , label = 'GI at S%s' %location[i])  # duration, mean_i
# show all storm events
plt.scatter(XX[cso_m,0], XX[cso_m,1], c='k',
                   marker=markers1[5])

plt.ylabel("Mean(i)(in/hr)")
plt.xlabel("Duration(hr)")
#plt.title("Duration vs Mean(i)")
plt.legend()
plt.grid()

plt.subplot(312)

plt.scatter(XX[CSO_it[0].astype(int),2], XX[CSO_it[0].astype(int),1], c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO')  # vol, mean_i
for i in range(n_sub):
    plt.scatter(XX[CSO_re[i],2], XX[CSO_re[i],1], c=color1[i],linewidths=0,
                   marker=markers1[0], s=s[i], alpha=Alpha, label = 'GI at S%s' %location[i])  # vol, mean_i
# show all storm events
plt.scatter(XX[cso_m,2], XX[cso_m,1], c='k',
                   marker=markers1[5])

plt.ylabel("Mean(i)(in/hr)")
plt.xlabel("Total Volume (in)")
plt.grid()

plt.subplot(313)
plt.scatter(XX[CSO_it[0].astype(int),3], XX[CSO_it[0].astype(int),1], c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO')  # max_i, mean_i
for i in range(n_sub): # show CSO events treated bgy GI
    plt.scatter(XX[CSO_re[i],3], XX[CSO_re[i],1], c=color1[i],linewidths=0,
                   marker=markers1[0], s=s[i], alpha=Alpha , label = 'GI at S%s' %location[i])  # max_i, mean_i
# show all storm events
plt.scatter(XX[cso_m,3], XX[cso_m,1], c='k',
                   marker=markers1[5])

plt.ylabel("Mean(i)(in/hr)")
plt.xlabel("Max(i)(in/hr)")
plt.grid()
plt.savefig("Fig2_%s_storm_clusters_vs_Mean_i_by_loc_cutoff_%s_DR_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

#%% Frequency Reduction - concentric circles
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 1
location = np.arange(1,6,1)
# XX[ length, i , volume, max_i]
cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
cso_m = CSO_it[0][cri1].astype(int)
cso_i_m = np.arange(len(CSO_it[0]))[cri1]  # the index of the storms satisfying the criteria
fig3= plt.figure(3, figsize=(8, 4))

ax = plt.gca()
ax.scatter(XX[CSO_it[0].astype(int),2]*DR/(bern_h), XX[CSO_it[0].astype(int),1]*tn*DR/(bern_h), c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO') 
ax.set_yscale('log')
ax.set_xscale('log')

for i in range(n_sub): # show CSO events treated by GI
    plt.scatter(XX[CSO_re[i],2]*DR/(bern_h), XX[CSO_re[i],1]*tn*DR/(bern_h), c=color1[i],linewidths=0, marker=markers1[0], s=s[i], alpha=Alpha ,
                label = 'GI at S%s' %location[i])  
plt.scatter(XX[cso_m,2]*DR/(bern_h), XX[cso_m,1]*tn*DR/(bern_h), c='k',
                   marker=markers1[5])
ax.set_yscale('log')
ax.set_xscale('log')

plt.xlabel("Tr")
plt.ylabel("Tn")
#
plt.xlim(0.1,10)
plt.ylim(.01,10)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="lower left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig3_%s_frequency reduction_by_loc_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

# Frequency reduction plot 2 - non-CSO, CSO, Loc-rlevent, Loc-irrelevant 


#%% Frequency Reduction - concentric circles tr/tn vs Tr
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 1
location = np.arange(1,6,1)
# XX[ length, i , volume, max_i]
cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
cso_m = CSO_it[0][cri1].astype(int)
cso_i_m = np.arange(len(CSO_it[0]))[cri1]  # the index of the storms satisfying the criteria
fig4= plt.figure(4, figsize=(8, 4))
ax = plt.gca()
ax.scatter(XX[CSO_it[0].astype(int),2]*DR/(bern_h), XX[CSO_it[0].astype(int),0]/tn, c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO') 
ax.set_yscale('log')
ax.set_xscale('log')

for i in range(n_sub): # show CSO events treated by GI
    plt.scatter(XX[CSO_re[i],2]*DR/(bern_h), XX[CSO_re[i],0]/tn, c=color1[i],linewidths=0, marker=markers1[0], s=s[i], alpha=Alpha ,
                label = 'GI at S%s' %location[i])  
plt.scatter(XX[cso_m,2]*DR/(bern_h), XX[cso_m,0]/tn, c='k',
                   marker=markers1[5])
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("tr/tn")
plt.xlabel("Tr")
plt.xlim(0.1,10)
plt.ylim(.1,100)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="lower right", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig4_%s_frequency reduction_tntr_Tr_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))
#%% 0220 - CSOn, nonCSO, Loc relevant, loc irrelevant 
# tr/tn (duration) vs Tr (mean intensity)
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 0.5
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# XX[ length, i , volume, max_i]
loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

cso_m = CSO_it[0].astype(int)

s1 = 40
# show non CSO storm events
fig5= plt.figure(5, figsize=(8, 4))
ax = plt.gca()
plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,0]/tn,s=s1, c='grey',
                   marker=markers1[0], alpha=0.2, edgecolor="grey", label = 'Non-CSO')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,0]/tn,s=s1, c='k',
                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'Large CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,0]/tn, c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,0]/tn,  c='r',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant')
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("tr/tn")
plt.xlabel("Tr")
plt.xlim(0.01,10)
plt.ylim(0.1,100)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig5_%s_CSO_Location_tntr_Tr_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))



#%% Frequency Reduction - concentric circles tr/tn vs Tn
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 1
location = np.arange(1,6,1)
# XX[ length, i , volume, max_i]
cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
cso_m = CSO_it[0][cri1].astype(int)
cso_i_m = np.arange(len(CSO_it[0]))[cri1]  # the index of the storms satisfying the criteria
fig6= plt.figure(6, figsize=(8, 4))
ax = plt.gca()
ax.scatter(XX[CSO_it[0].astype(int),1]*tn*DR/(bern_h), XX[CSO_it[0].astype(int),0]/tn, c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO') 
ax.set_yscale('log')
ax.set_xscale('log')

for i in range(n_sub): # show CSO events treated by GI
    plt.scatter(XX[CSO_re[i],1]*tn*DR/(bern_h), XX[CSO_re[i],0]/tn, c=color1[i],linewidths=0, marker=markers1[0], s=s[i], alpha=Alpha ,
                label = 'GI at S%s' %location[i])  
plt.scatter(XX[cso_m,1]*tn*DR/(bern_h), XX[cso_m,0]/tn, c='k',
                   marker=markers1[5])
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("tr/tn")
plt.xlabel("Tn")
plt.xlim(0.1,10)
plt.ylim(0.1,100)
plt.rcParams.update({'font.size': 12})
plt.legend()
#plt.legend(loc="lower left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig6_%s_frequency reduction_tntr_Tn_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

#%% 0220 - CSOn, nonCSO, Loc relevant, loc irrelevant 
# tr/tn (duration) vs Tn (mean intensity)
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 0.5
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# XX[ length, i , volume, max_i]
loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

cso_m = CSO_it[0].astype(int)

s1 = 40
# show non CSO storm events
fig7= plt.figure(7, figsize=(8, 4))
ax = plt.gca()
plt.scatter(XX[nonCSO,1]*tn*DR/(bern_h), XX[nonCSO,0]/tn,s=s1, c='grey',
                   marker=markers1[0], alpha=0.2, edgecolor="grey", label = 'Non-CSO')
# show CSO storms
plt.scatter(XX[CSOs,1]*tn*DR/(bern_h), XX[CSOs,0]/tn,s=s1, c='k',
                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,1]*tn*DR/(bern_h), XX[loc_irr,0]/tn, c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,1]*tn*DR/(bern_h), XX[loc_rel,0]/tn,  c='r',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="g",
                label = 'Loc-relevant CSOs')
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("tr/tn")
plt.xlabel("Tn")
plt.xlim(0.01,10)
plt.ylim(0.1,100)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="upper right", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig7_%s_CSO_Location_tntr_Tn_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

#%% 0220 - CSOn, nonCSO, Loc relevant, loc irrelevant 
# tr/tn (duration) vs Tr (mean intensity)
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 0.5
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# XX[ length, i , volume, max_i]
loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

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
fig10= plt.figure(10, figsize=(8, 6))
ax = plt.gca()
plt.plot([GN,10],GN*np.ones(2), c='b',label = 'GN-Lower Bound')
plt.plot(GN*np.ones(2),[GN,10], c ='b')
plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='grey',
                   marker=markers1[0], alpha=0.2, edgecolor="grey", label = 'Non-CSOs')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'Persisting CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant CSOs')
plt.plot(Tr,env_y,label = 'GN-Envelop',c='g')
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("Tn")
plt.xlabel("Tr")
plt.xlim(0.01,10)
plt.ylim(0.01,10)
plt.rcParams.update({'font.size': 16})
#plt.legend()
plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig10_GN_%s_tn_%s_Tn_Tr_dr_%s.png" %(city,tn,DR))

#%% Adjust the envelope lines  
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 0.5
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
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
a = 0.25 # adjustion fraction for the storage increase from infiltration of long duration storms
# only appied to the horozontal bound 
env_y = []
Gn_x=GN
Gn_y=GN*a
Gn=GN
Tr = np.arange(Gn+0.001,6,0.01)

for Tr1 in Tr:
    zz = -(np.exp(-Tr1/Gn)*Tr1/Gn)
#    print(zz)
    w = lambertw(zz,0)
    Tn1 = (Gn*Tr1)/(Tr1+Gn*w.real)
    env_y.append(Tn1-GN*(1-a))
#plt.plot(env_y)
# show non CSO storm events
fig11= plt.figure(11, figsize=(8, 6))
ax = plt.gca()
plt.plot([Gn_x,10],(Gn_y)*np.ones(2), c='b',label = 'GN-Adjusted Lower Bound')   # y = HN
plt.plot(Gn_x*np.ones(2),[Gn_y,10], c ='b')                    # x = HN
plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='grey',
                   marker=markers1[0], alpha=0.2, edgecolor="grey", label = 'Non-CSOs')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'Persisting CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant CSOs')
plt.plot(Tr,env_y,label = 'GN-Ajusted Envelop',c='g')
ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("Tn")
plt.xlabel("Tr")
plt.xlim(0.01,10)
plt.ylim(0.01,10)
plt.rcParams.update({'font.size': 16})
#plt.legend()
#plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()

plt.grid()
plt.savefig("Fig11_GN_%s_tn_%s_Tn_Tr_dr_%s_adjusted_envelope_line_storage_adjusted.png" %(city,tn,DR))

#%% Adjust the envelope linesby multiplying a fraction
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 0.5
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
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
a = 0.25 # adjustion fraction for the storage increase from infiltration of long duration storms
# only appied to the horozontal bound 
env_y = []
Gn_x=GN
Gn_y=a*GN
Gn=GN
Tr = np.arange(Gn+0.0001,10,0.01)

for Tr1 in Tr:
    zz = -(np.exp(-Tr1/Gn)*Tr1/Gn)
#    print(zz)
    w = lambertw(zz,0)
    Tn1 = a*(Gn*Tr1)/(Tr1+Gn*w.real)
    env_y.append(Tn1)
#plt.plot(env_y)
# show non CSO storm events
fig13= plt.figure(13, figsize=(8, 6))
ax = plt.gca()
#plt.plot([Gn_x,10],(Gn_y)*np.ones(2), c='b',label = 'GN-Adjusted Lower Bound')   # y = HN
#plt.plot(Gn_x*np.ones(2),[Gn_y,10], c ='b')                    # x = HN
#plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='grey',
#                   marker=markers1[0], alpha=0.2, edgecolor="grey", label = 'Non-CSOs')
## show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'Persisting CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='b',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant CSOs')
plt.plot(Tr,env_y,label = 'GN-Ajusted Envelop',c='g')

#
#qcso= cutoff/(1-DR*gi_s_p) # with GI, the Qcso threshold is adjusted based on the area treated
#HN = tn*qcso*DR/bern_h
#env_y = []
#Hn_x=HN
#Hn_y=HN*a
#Hn=HN
#Tr = np.arange(Hn+0.0001,10,0.01)
#for Tr1 in Tr:
#    zz = -(np.exp(-Tr1/Hn)*Tr1/Hn)
##    print(zz)
#    w = lambertw(zz,0)
#    Tn1 = (Hn*Tr1)/(Tr1+Hn*w.real)
#    env_y.append(Tn1-HN*(1-a))
#    
#plt.plot(Tr,env_y,label = 'GN-Envelop_w/GI',c='g',  linestyle='-.')


ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("Tn")
plt.xlabel("Tr")
plt.xlim(0.01,10)
plt.ylim(0.01,10)
plt.rcParams.update({'font.size': 16})
#plt.legend()
#plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()

plt.grid()
plt.savefig("Fig12_GN_%s_tn_%s_Tn_Tr_dr_%s_HN_line_adjusted.png" %(city,tn,DR))
#
print'total cso no. = ', len(CSO_it[0])
print 'loc-irr = ', len(loc_irr)
print 'loc-rel = ', len(loc_rel)
print 'CSOs = ', len(CSOs)
print 'non-CSOs', len(nonCSO)