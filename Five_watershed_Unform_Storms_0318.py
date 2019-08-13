"""
Created on Sun May 21 20:42:17 2017

@author: Fengwei Hung
"""
from __future__ import division
import os
os.chdir('C:/School/2nd Paper/model')  # for my computer only
import numpy as np
import matplotlib.pyplot as plt
from cat_runoff_w_infil import sub_runoff
from cat_runoff_w_infil import Muskingun
import pandas as pd
from def_storm_CDO import storm_id
from rainfall_analysis_hour import storms
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
bern_h = 6 #in

###########################
##  Indepent experiment
###########################
### Chaning parameters 
location = 0
L = 100            # length of the watershed
city = 'Philadelphia'

P_crit = 0.001 # in/hr
P_lag = 2
metric = 'in' # indicate the metric in the input file: mm or in
qcso= cutoff
HN = tn*qcso*DR/bern_h
#%% rainfall generater
r_b = 40
rain = []
for v in np.arange(0.4,2,0.1):
#for i in np.arange(0.2,0.41,0.01):
    for tr in np.arange(1,12,1):
        i = v/tr
        rainfall = np.zeros(int((tr+r_b)))
        rainfall[0:int(tr)] = i*np.ones(int(tr))  #storm duration 12 hour
        rain=rain+list(rainfall)
        
rng = pd.date_range('1/1/2000', periods=len(rain), freq='H')   
end_time = rng[len(rain)-1]
print('End Data = %s' %end_time)
data = pd.Series(rain, index=rng)   
df = pd.DataFrame(rain, index=rng, columns=['P'])
#writer = pd.ExcelWriter('%s_uniform_stroms.xlsx' %city)
#df.to_excel(writer,'P')
filename = '%s_uniform_stroms.csv' %city  # unit: =in/hr
df.to_csv(filename, sep=',', encoding='utf-8')
  
run_dates =     ['2000-01-01', '2001-08-29']
#run_dates =     ['2000-01-01', '2001-06-24']

n_yr = 1
df_P1, data_s = storm_id(filename, run_dates,P_crit, P_lag,metric,cutoff)   # mean_i, length, max_i, storm_id, volume(mm)
df_P = df_P1.resample('10min').bfill()
rain= df_P['P'].values.tolist()
rng2 = df_P.index.values
plt.plot(rain)
##############
#### plot ####
##############
#%% 
color1 =['b','g','r','m','y','c','k','b','g','r']
line=['-.','--',':','-.']
### parameters for Muskingun

cso_id =np.ones(len(data_s['i']))*-1
data_s['cso_id'] = cso_id
CSO_it = []  ## index of the storms that generate CSO
CSO_re = []  ## CSO event reduction
CSO_t=[]
Peak_t=[]
CSO_f=[]
CSO_V = []
err_v = 0.001  # the monimum flow exceeding q_CSO (in) for a CSO event/ used for aoiding numberical errors
for location in range (0,6):
    k = tn/n_sub
    x = 0.10 
    c0=(-k*x+0.5*dt)/(k-k*x+0.5*dt)
    c1=(k*x+0.5*dt)/(k-k*x+0.5*dt)
    c2=(k-k*x-0.5*dt)/(k-k*x+0.5*dt)
    End_T= len(rain)*dt
    Q = np.zeros(len(rain))
    Q_all = []    
    x1= np.arange(0,End_T+dt,dt)
    Ft = 0
    for j in range(1,(n_sub+1)):  # calculate hydrograph at the each subcatchment
        if location == j:
            gi_s = all_area*gi_s_p     # Surface area of GI (ft^2)   
            gi_d = gi_s*DR            # drainage area (ft^2)  
            w = sub_area[j-1]/L          # width (ft)
            w_gi = gi_d/L                # width
            rf_s1, f_s1, d_s1 =  sub_runoff(rain,sub_area[j-1],imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h) # [ p, imp, gi]
            rf = rf_s1[0]+ rf_s1[1]+ rf_s1[2]
            Q = Muskingun(Q, rf,dt,End_T,c0,c1,c2,rain)            
            Q_t=np.array(Q)
            Q_t[Q_t<0]=0
            Q_all.append(Q_t)
            F =np.sum(np.array(f_s1[1])*dt)/12*gi_s*7.48/1000  # kilo gallon       
        else:
            gi_s = 0     # Surface area of GI (ft^2)   
            gi_d = 0     # drainage area (ft^2)  
            w = sub_area[j-1]/L          # width (ft)
            w_gi = gi_d/L                # width
            rf_s1, f_s1, d_s1 = sub_runoff(rain,sub_area[j-1],imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h) # [ p, imp, gi]
            rf = rf_s1[0]+ rf_s1[1]+ rf_s1[2]
            Q = Muskingun(Q, rf,dt,End_T,c0,c1,c2,rain)    
            Q_t=np.array(Q)
            Q_t[Q_t<0]=0
            Q_all.append(Q_t)  ## Q from S1 to S5
            F =np.sum(np.array(f_s1[1])*dt)/12*gi_s*7.48/1000  # kilo gallon
        Ft += F
        print "Subcatchment = %s" %j
#    plt.figure(figsize=(8, 6))
##    plt.plot(Q_t)
#    plt.plot(Q_t[1120:1250])
#    plt.grid()
#    plt.tight_layout()
#    plt.savefig("Fig11_uniform_storms_discharges_loc_%s.png" %location)
        
    #### CSO Estimation 

    overflow=0  # the amount of overflow (cfs)

    cso_i=[]     # the storm_id that indicates which storm generates the CSO
    q_fl=0  # the flow causes CSO 
    Q_fl= [] # the amount of overflow for each storm (cf)
    storm_i = 0     #
    cso_vol = 0     # the total CSO volume of a storm event           

    for ii in range(len(Q_t)):
        if df_P['storm_id'][ii]!=0:
            storm_i = df_P['storm_id'][ii]            
                       
        if (Q[ii]>(cutoff+err_v)):
            cso_i.append(storm_i)
            overflow= (Q_t[ii]-cutoff)*dt*3600
            Q_fl.append(overflow)  
#            print('Q = ', Q[ii])
        else:
            cso_i.append(0)
            Q_fl.append(0)  
            
    df_P['cso_i'] = cso_i  

    
    Q_fl_a = np.array(Q_fl)             # array of overflow (volume)
#    Q_fl_i = np.where(Q_fl_a>0)         # the index of the overflow
#    cso_b = np.array(cso_i)[Q_fl_i]     # find the storm index of the CSOs
#    cso_b = np.append(cso_b,[0])        # add 0
#    Q_fl_b = Q_fl_a[Q_fl_i]             # overflows only
    cso_a = []
    cso_v = []
# calcluate CSO volume for individual storms    
    
    s_id = 0
    for i in np.arange(0,len(Q_fl)-1):
        if df_P.loc[rng2[i],'storm_id']>0:
            s_id = df_P.loc[rng2[i],'storm_id']
            cso_vol = cso_vol + Q_fl_a[i]*7.48/(1000) # kgal
            if df_P.loc[rng2[i+1],'storm_id']==0:
                df_P.loc[rng2[i+1],'storm_id'] = s_id
                if i == len(Q_fl)-2:
                    cso_v.append(cso_vol)
                    if cso_vol>0.0001:
                        cso_a.append(int(s_id))
                    print('s_id = %s' %s_id)
                    print('cso_vol = %s' %cso_vol)
                    cso_vol=0                    
            elif df_P.loc[rng2[i],'storm_id'] < df_P.loc[rng2[i+1],'storm_id'] or i == len(Q_fl)-2:
                cso_v.append(cso_vol)
                if cso_vol>0.0001:
                    cso_a.append(int(s_id))
#                print('s_id = %s' %s_id)
#                print('cso_vol = %s' %cso_vol)
                cso_vol=0
    CSO_V.append(np.array(cso_v))
    cso_f=len(cso_a)     ## count number of CSO event
    CSO_f.append(cso_f)   
    cso_ii = np.array(cso_a)-1     ## this is the index for assigning location the cso_id column 
    cso_id[cso_ii] = location
    CSO_it.append(cso_ii)       ## CSO_it are indeces (not the storm ids), too
    if location != 0:
        CSO_re.append(set(CSO_it[0])-set(CSO_it[location]))

    else:                       # cso_id indicate cso events in data_S
#        df_P['Q']=Q
        for j in range(len(cso_ii)):
            data_s['cso_id'][cso_ii[j]] = 1
          
        print "No. of CSO events: %s" %cso_f
    cso_kgal=sum(np.array(cso_v))  # k-gal
    CSO_t.append(cso_kgal)
    print "cso = %s kgal" %cso_kgal
    peak = max(Q_t)
    Peak_t.append(peak)
    print "Peak flow = %s" %peak
        
#data_cso = np.zeros(len(cso_id))
data_s['CSO_V'] = CSO_V[0]
#data_cso_re = np.zeros(len(cso_id))
#data_cso[CSO_it[0].astype(int)] = CSO_V[0]

    #%%  save to file
#import pickle
writer = pd.ExcelWriter('%s_Uniform_Storm_dataframe.xlsx' %city)
data_s.to_excel(writer,'Data_s')

# CSO_t
d_cso = np.zeros([6,2])
d_cso[0:6,0]=CSO_t
d_cso[0:6,1]=CSO_f

data_CSO =pd.DataFrame(data = d_cso, columns = ["CSO_Kgal","CSO_f"])
data_CSO.to_excel(writer,'Annual_CSO')

for i in range(n_sub+1):    
    data_CSO_it = pd.DataFrame(data =CSO_it[i], columns = ["CSO_%s" %i])
    data_CSO_it.to_excel(writer, "CSO_it_%s" %i)

for i in range(n_sub+1):    
    data_CSO_V = pd.DataFrame(data =CSO_V[i], columns = ["CSO_%s" %i])
    data_CSO_V.to_excel(writer, "CSO_V%s" %i)
#CSO_it.to_excel(writer, 'Sheet2')
writer.save()


   
#%%    
markers1 = ['o','>','<','1','x','.','1','*','|','_',' ']
s = np.array([500,350,250,160,80,20])
location = np.arange(1,6,1)

N = len(data_s['i'])
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s['i'][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# convert CSO_re (set) to lists

CSO_V_re=[]
for i in range(n_sub):
    k=0
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
#    plt.xlim(0.6,10)
#    plt.ylim(0.1,10)
    plt.rcParams.update({'font.size': 14})
#    plt.legend(loc="lower left", markerscale=1., scatterpoints=1, fontsize=10)   
#    lgnd = plt.legend(loc="upper left", scatterpoints=1, fontsize=13)
#    lgnd.legendHandles[0]._sizes = [50]
#    lgnd.legendHandles[1]._sizes = [50]
#    lgnd.legendHandles[2]._sizes = [50]
#    lgnd.legendHandles[3]._sizes = [50]   

 #change the marker size manually for both lines

    plt.tight_layout()
    plt.grid()
    plt.savefig("Fig%s_%s_TrTn_S_%s.png" %((location[loc],city,location[loc])))





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

plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(12, 12))
location = np.arange(1,6,1)
Alpha =  1 # transparebcy

plt.subplot(311)

#plt.figure(figsize=(8, 4))  # Tr : Tn
plt.scatter(XX[CSO_it[0].astype(int),0], XX[CSO_it[0].astype(int),1], c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO')  # vol, mean_i
for i in range(n_sub):
    plt.scatter(XX[CSO_re[i],0], XX[CSO_re[i],1], c=color1[i],linewidths=0,
                   marker=markers1[0], s=s[i], alpha=Alpha , label = 'GI at S%s' %location[i])  # duration, mean_i
# show all storm events
plt.scatter(XX[cso_m,0], XX[cso_m,1], c='k',
                   marker=markers1[5])
#plt.xlim(0.3,10)
#plt.ylim(0.15,0.4)
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
plt.savefig("Fig6_%s_storm_clusters_vs_Mean_i_by_loc_cutoff_%s_DR_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

#%% Frequency Reduction - concentric circles
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 1
location = np.arange(1,6,1)
# XX[ length, i , volume, max_i]
cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
cso_m = CSO_it[0][cri1].astype(int)
cso_i_m = np.arange(len(CSO_it[0]))[cri1]  # the index of the storms satisfying the criteria
plt.figure(figsize=(8, 4))  # Volume vs Maximum Intensity
ax = plt.gca()
ax.scatter(XX[CSO_it[0].astype(int),2]*DR/(bern_h), XX[CSO_it[0].astype(int),1]*tn*DR/(bern_h), c='grey',linewidths=0,
               marker=markers1[0], s=s[4], alpha=Alpha , label = 'CSO') 
#ax.set_yscale('log')
#ax.set_xscale('log')

for i in range(n_sub): # show CSO events treated by GI
    plt.scatter(XX[CSO_re[i],2]*DR/(bern_h), XX[CSO_re[i],1]*tn*DR/(bern_h), c=color1[i],linewidths=0, marker=markers1[0], s=s[i], alpha=Alpha ,
                label = 'GI at S%s' %location[i])  
plt.scatter(XX[cso_m,2]*DR/(bern_h), XX[cso_m,1]*tn*DR/(bern_h), c='k',
                   marker=markers1[5])

plt.grid()
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.xlabel("Tr")
plt.ylabel("Tn")
#
#plt.xlim(1,2)
#plt.ylim(1,5)
plt.rcParams.update({'font.size': 12})
#plt.legend()
#plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()

plt.savefig("Fig7_%s_frequency reduction_by_loc_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

# Frequency reduction plot 2 - non-CSO, CSO, Loc-rlevent, Loc-irrelevant 


#%% Frequency Reduction - concentric circles tr/tn vs Tr
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 1
location = np.arange(1,6,1)
# XX[ length, i , volume, max_i]
cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
cso_m = CSO_it[0][cri1].astype(int)
cso_i_m = np.arange(len(CSO_it[0]))[cri1]  # the index of the storms satisfying the criteria
plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
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
#plt.xlim(0.3,10)
#plt.ylim(.1,10)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="lower right", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig8_%s_frequency reduction_tntr_Tr_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))
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
plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
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
#plt.xlim(0.4,10)
#plt.ylim(0.1,10)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig8_%s_CSO_Location_tntr_Tr_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))



#%% Frequency Reduction - concentric circles tr/tn vs Tn
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 1
location = np.arange(1,6,1)
# XX[ length, i , volume, max_i]
cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
cso_m = CSO_it[0][cri1].astype(int)
cso_i_m = np.arange(len(CSO_it[0]))[cri1]  # the index of the storms satisfying the criteria
plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
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
#plt.xlim(0.1,10)
#plt.ylim(0.1,10)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="lower left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig9_%s_frequency reduction_tntr_Tn_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

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
plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
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
#plt.xlim(0.01,0.3)
#plt.ylim(0.1,80)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="lower left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig9_%s_CSO_Location_tntr_Tn_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

#%% 0220 - CSOn, nonCSO, Loc relevant, loc irrelevant 
# tr/tn (duration) vs Tr (mean intensity)
from scipy.special import lambertw

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
#Tn = np.arange(1,4,0.1)
Tr = np.arange(HN+0.01,6,0.01)
env_y = []
Hn=HN
for Tr1 in Tr:
    zz = -(np.exp(-Tr1/Hn)*Tr1/Hn)
#    print(zz)
    w = lambertw(zz,0)
    Tn1 = (Hn*Tr1)/(Tr1+Hn*w.real)
    env_y.append(Tn1)


plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
ax = plt.gca()
plt.plot([HN,10],HN*np.ones(2), c='b',label = 'HN-Line')
plt.plot(HN*np.ones(2),[HN,10], c ='b')
plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='grey',
                   marker=markers1[0], alpha=0.2, edgecolor="grey", label = 'Non-CSO')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant CSOs')
plt.plot(Tr,env_y,label = 'HN-Envelop')

ax.set_yscale('log')
ax.set_xscale('log')

plt.ylabel("Tn")
plt.xlabel("Tr")
plt.xlim(0.1,4)
plt.ylim(0.1,4)
plt.rcParams.update({'font.size': 12})
plt.legend()
plt.legend(loc="upper left", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig10_%s_tn_%s_Tn_Tr_dr_%s.png" %(city,tn,DR))

#%% 0220 - CSOn, nonCSO, Loc relevant, loc irrelevant 
#Tr = np.arange(1.21,6,0.01)
#env_y = []
#
#for Tr1 in Tr:
#    zz = -(np.exp(-Tr1/HN)*Tr1/HN)
##    print(zz)
#    w = lambertw(zz,0)
#    Tn1 = (HN*Tr1)/(Tr1+HN*w.real)
#    env_y.append(Tn1)


plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
ax = plt.gca()
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.plot([HN,10],HN*np.ones(2), c='b',label = 'HN-Line')
plt.plot(HN*np.ones(2),[HN,12], c='b')
# show CSO storms
plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'CSOs')
# plot Loc-irrelevent storms
plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
                label = 'Loc-irrelevant CSOs')  
plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="r",
                label = 'Loc-relevant CSOs')
plt.plot(Tr,env_y,label = 'HN-Envolope')

    

plt.rcParams.update({'font.size': 12})
plt.legend()
plt.legend(loc="upper left", scatterpoints=1, fontsize=12)
plt.xlim(0,2)
plt.ylim(0,4)
plt.legend()
plt.legend(loc="upper right", scatterpoints=1, fontsize=12)
plt.tight_layout()

plt.grid()
plt.savefig("Fig11_%s_tn_%s_Tn_Tr_dr_%s.png" %(city,tn,DR))
