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
from scipy.special import lambertw

#from Hydrograph_plot import Muskingun_plot 
cutoff = 0.10 # Q_CSO: cfs
i_metric = '90perc_i'

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
###########################
##  Indepent experiment
###########################
### Chaning parameters 
location = 0
L = 100            # length of the watershed
#city = 'Philadelphia'
city = 'Seattle'
filename = city +'_1980_2013.csv'  # unit: =in/hr

P_crit = 0.01 # in/hr
P_lag = 2
metric = 'in' # indicate the metric in the input file: mm or in

#%%
End_year = 1981
run_dates =     ['1980-01-01', str(End_year)+'-12-31']
n_yr = End_year-1980
df_P1, data_s = storm_id(filename, run_dates,P_crit, P_lag,metric)   # mean_i, length, max_i, storm_id, volume(mm)
df_P = df_P1.resample('10min').bfill()
rain= df_P['P'].values.tolist()


#%% 
##############
#### plot ####
##############
color1 =['b','g','r','m','y','c','k','b','g','r']
line=['-.','--',':','-.']
### parameters for Muskingun
tn = 3

cso_id =np.ones(len(data_s['storm_id']))*-1
data_s['cso_id'] = cso_id
CSO_it = []  ## index of the storms that generate CSO
CSO_re = []  ## CSO event reduction
CSO_t=[]
Peak_t=[]
CSO_f=[]
CSO_V = []
err_v = 0.001  # the monimum flow exceeding q_CSO (in) for a CSO event/ used for aoiding numberical errors
F_V = []
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
            rf = rf_s1[0]+ rf_s1[1]+ rf_s1[2]  #ft3
            ft = f_s1[0]+f_s1[1]               #in/hr 
            Q = Muskingun(Q, rf,dt,End_T,c0,c1,c2,rain)            
            Q_t=np.array(Q)
            Q_t[Q_t<0]=0
            Q_all.append(Q_t)
            F =np.array(ft)  # in       
        else:
            gi_s = 0     # Surface area of GI (ft^2)   
            gi_d = 0     # drainage area (ft^2)  
            w = sub_area[j-1]/L          # width (ft)
            w_gi = gi_d/L                # width
            rf_s1, f_s1, d_s1 = sub_runoff(rain,sub_area[j-1],imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h) # [ p, imp, gi]
            rf = rf_s1[0]+ rf_s1[1]+ rf_s1[2]
            ft = f_s1[0]+f_s1[1]

            Q = Muskingun(Q, rf,dt,End_T,c0,c1,c2,rain)    
            Q_t=np.array(Q)
            Q_t[Q_t<0]=0
            Q_all.append(Q_t)  ## Q from S1 to S5
            F =np.array(ft)         # in
        Ft += F*dt
        print "Subcatchment = %s" %j
        
    #### CSO Estimation 

    overflow=0  # the amount of overflow (cfs)

    cso_i=[]     # the storm_id that indicates which storm generates the CSO
    
    q_fl=0  # the flow causes CSO 
    Q_fl= [] # the amount of overflow for each storm (cfs)
    storm_i = 0     #
    cso_vol = 0     # the total CSO volume of a storm event     
    f_vol = 0       # the total infiltration (kgal)      

    cso_i_only = [] # the storm id that generates CSO (excluded 0)
    Q_fl_only =[]   # the amount of overflow for each storm (cfs; excluded 0)
    Ft_only = []    # the amount of infiltration during CSO events
    for ii in range(len(Q_t)):
        if df_P['storm_id'][ii]!=0:
            storm_i = df_P['storm_id'][ii]            
                       
            if (Q[ii]>(cutoff+err_v)):
                cso_i.append(storm_i)
                overflow= (Q_t[ii]-cutoff)*dt*3600
                Q_fl.append(overflow)  # Q_overflow
    #            print('Q = ', Q[ii])
                cso_i_only.append(storm_i)
                Q_fl_only.append(overflow)
                Ft_only.append(Ft[ii])
            else:
                cso_i.append(0)
                Q_fl.append(0)  
            
#    df_P['cso_i'] = cso_i  
#    df_P['Infil'] = Ft  
    Q_fl_a = np.array(Q_fl)             # array of overflow (volume)
    cso_a = []
    cso_v = []
    f_v = []            # the infiltrated water volume (kgal) for all storms
# calcluate CSO volume for individual storms   
    cso_i_only.append(0)
    cso_i_only_set = set(cso_i_only)

    
    for i in np.arange(0,len(cso_i_only)-1):
        if cso_i_only[i] == cso_i_only[i+1]:     # cso vol for an event
            cso_vol = cso_vol + Q_fl_only[i]*7.48/(1000) # kgal
            f_vol = f_vol + Ft_only[i]                      #kgal
        else:  
            cso_a.append(int(cso_i_only[i]))                     # save the storm id; transition to a new storm
            cso_vol = cso_vol + Q_fl_only[i]*7.48/(1000) # kgal
            cso_v.append(cso_vol)
            f_v.append(f_vol)
            cso_vol = 0
            f_vol = 0    
#    s_id = 0
#    for i in np.arange(0,len(Q_fl)-1):
#        if df_P['storm_id'][i]>0:
#            s_id = df_P['storm_id'][i]
#            cso_vol = cso_vol + Q_fl_a[i]*7.48/(1000) # kgal
#            f_vol = f_vol + Ft[i]                      #kgal
#            print('s_id = %s' %s_id)
#            print('cso_vol = %s' %cso_vol)
#            if df_P['storm_id'][i+1]==0:   
#                df_P['storm_id'][i+1] = s_id          # a new storm
#                if i == len(Q_fl)-2:                  # the end of the data
#                    cso_v.append(cso_vol)
#                    f_v.append(f_vol)
#                    if cso_vol>0.0001:
#                        cso_a.append(int(s_id))
#
#                cso_vol=0      
#                f_vol = 0
#            elif df_P['storm_id'][i] < df_P['storm_id'][i+1] or i == len(Q_fl)-2:
#                cso_v.append(cso_vol)
#                f_v.append(f_vol)
#                if cso_vol>0.0001:
#                    cso_a.append(int(s_id))
##                print('s_id = %s' %s_id)
##                print('cso_vol = %s' %cso_vol)
#                cso_vol=0
#                f_vol = 0

    CSO_V.append(np.array(cso_v))
    F_V.append(f_v) 
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
          
#        print "No. of CSO events: %s" %cso_f
    cso_kgal=sum(np.array(cso_v))  # k-gal
    CSO_t.append(cso_kgal)
#    print "cso = %s kgal" %cso_kgal
    peak = max(Q_t)
    Peak_t.append(peak)
#    print "Peak flow = %s" %peak
        
#data_cso = np.zeros(len(cso_id))   
data_s['CSO_V'] = np.ones(len(data_s['storm_id']))*0
data_s['Infil'] = np.ones(len(data_s['storm_id']))*0
data_s.iloc[CSO_it[0],8] = CSO_V[0] # CSO volume without GI
data_s.iloc[CSO_it[4],9] = F_V[4]   # infiltration of the persisting CSOs

#data_cso_re = np.zeros(len(cso_id))
#data_cso[CSO_it[0].astype(int)] = CSO_V[0]

    #%%  save to file
writer = pd.ExcelWriter('%s_Storm_dataframe_i_adjusted.xlsx' %city)
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

N = len(data_s['storm_id'])
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s[i_metric][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# convert CSO_re (set) to lists
count = 0
RE_i = []  # nested list of index for CSO reduction calculation
for i in range(n_sub):
    re_i = []  # index for calculating CSO reduction
    count = 0
    for j in CSO_it[0]:
        if j in CSO_it[i+1]:
            re_i.append(count)
        count += 1
    RE_i.append(re_i)   
    


CSO_V_re=[]  # 
for i in range(n_sub):
    CSO_re[i] = map(int, list(CSO_re[i]))  # CSO stroms that were fully treated
    cso_re = np.array(data_s.iloc[CSO_it[0],8])
    count_1 = 0  # the index of all CSO storms (the case of no GI)
#    count_2 = 0  # the index of the persisting CSOs storms
    for j in RE_i[i]:
        cso_re[j] = cso_re[j] - CSO_V[i+1][count_1]
        count_1 += 1
    CSO_V_re.append(cso_re)   
    


    
#%% plot no 3: focus on manageable storms
color1 =['b','y','r','y','c','k','m','g','r','lightgreen']

Alpha = 0.8
S_gi = gi_s_p*bern_h/12*7.48*43560/1000 # 5.43 kgal GI's storage
# remove extreme storms

## GI at S1   
location=np.arange(1,6,1)
#cri1 = (XX[CSO_it[0].astype(int),2]<100)*(XX[CSO_it[0].astype(int),2]>0.0)
#cso_m = CSO_it[0].astype(int)
#cso_i_m = np.arange(len(CSO_it[0]))  # the index of the storms satisfying the criteria


for loc in np.arange(0,n_sub,1):

    CSO_V_RE = CSO_V_re[loc]
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
                   marker=markers1[0], s=CSO_V[0]*10, alpha=Alpha , label = 'CSO Volume') 

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
    plt.ylim(0.1,10)
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
#    for i in np.arange(0,9):   # loop hthough years        
        storm_iy=[]                # index of storms for a particular year i
        for j in CSO_it[location].astype(int):  # for CSO stroms j
            if years[i] == data_s.index.year[j]:  #the year of the storm j = year i
                storm_iy.append(j)   # stored storm j in storm_iy
#                print(j)

        storm_ay.append(storm_iy)        #identify the storm # in scenarios that GI sitted at each subwatershed
        
# calculate annual CSO volume        
        ns = len(storm_ay[i])       # the number of CSO storms of an year
#        print(ns)
        cso_y[i] = np.array(CSO_V[location][storm_iy]).sum()  # kgal
#        print(CSO_V[location][n0:(n0+ns)])
        cso_n[i] = ns
        n0=n0+ns
#        print('n0= ', n0)
#        print('Year =', 1980+i)
#        print('CSO volume',cso_y[i])
#        print('storm_iy=',storm_iy)
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
#np.corrcoef(cso_y_l[1],cso_y_l[3]) 
#np.corrcoef(cso_y_l[2],cso_y_l[3]) 
#np.corrcoef(cso_y_l[3],cso_y_l[4]) 
#np.corrcoef(cso_y_l[4],cso_y_l[5]) 

#%% Read csv files 
v_temp = pd.read_csv('cso_y_l.csv', sep=',',header=None)
cso_y_l = v_temp.values
n_temp = pd.read_csv('cso_y_n.csv', sep=',',header=None)
cso_y_n = n_temp.values
#%% Plot CSO reduction

CSO_yrm=[]
for i in np.arange(1,6):
    cso_y1 =cso_y_l[0]-cso_y_l[i]
    CSO_yrm.append(cso_y1)
    
CSO_reV = np.average(CSO_yrm, axis = 1)  # average annual CSO volume
CSO_V_se = 1.96*np.std(CSO_yrm,axis =1)/(n_yr**0.5) # standard error of annual CSO volume

CSO_yrn = []
for i in np.arange(1,6):
    cso_y1 =cso_y_n[0]-cso_y_n[i]
    CSO_yrn.append(cso_y1)
CSO_reN = np.average(CSO_yrn, axis = 1)  # average annual CSO volume    
CSO_N_se = 1.96*np.std(CSO_yrn,axis =1)/(n_yr**0.5) # standard error of annual CSO volume
    
fig1 = plt.rcParams.update({'font.size': 14})
width = 0.5  
fig4= plt.figure(4, figsize=(6, 3))
plt.bar(range(5), CSO_reV, width, color='b', yerr=CSO_V_se,alpha =0.5)
plt.xlim(-0.2,5.)
plt.ylim(0,30)
plt.ylabel("Volume Reduction (kgal/ac)")
plt.xlabel("GI Location")
plt.xticks( np.arange(0.,4.2,1), ('S1', 'S2', 'S3', 'S4', 'S5'), rotation='vertical' )

plt.grid()
plt.tight_layout()
plt.savefig("Fig1_%s_CSO_vol_re_dr_%s_bh_%s.png" %(city,DR,bern_h))

width = 0.5  
fig5= plt.figure(5, figsize=(6, 3))
#plt.plot(range(7), np.repeat(CSO_yn[0],7))
plt.bar(range(5), CSO_reN, width, color='b', yerr=CSO_N_se,alpha =0.5)
plt.xlim(-0.2,5)
plt.ylim(0,8)

plt.ylabel("Occurrance Reduction")
plt.xlabel("GI Location")
#plt.yticks( np.arange(min(CSO_f)-1, max(CSO_f)+1, 1))
plt.xticks(  np.arange(0.,4.2,1), ('S1', 'S2', 'S3', 'S4', 'S5'), rotation='vertical' )
#plt.xticks(x, labels, rotation='vertical')
plt.grid()
plt.tight_layout()
plt.savefig("Fig1_%s_CSO_occurance_re_dr_%s_bh_%s.png" %(city,DR,bern_h))



#%% 0220 - CSOn, nonCSO, Loc relevant, loc irrelevant 
# tr/tn (duration) vs Tn (mean intensity)
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 0.5
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s[i_metric][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# XX[ length, i , volume, max_i]
loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

cso_m = CSO_it[0].astype(int)

s1 = 40
# show non CSO storm events
fig9= plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
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
plt.ylim(0.1,80)
plt.rcParams.update({'font.size': 12})
#plt.legend()
plt.legend(loc="upper right", scatterpoints=1, fontsize=12)

plt.tight_layout()


plt.grid()
plt.savefig("Fig9_%s_CSO_Location_tntr_Tn_cutoff_%s_dr_%s_bh_%s.png" %(city,cutoff,DR,bern_h))

#%% 0220 - CSOn, nonCSO, Loc relevant, loc irrelevant 
# tr/tn (duration) vs Tr (mean intensity)
color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
Alpha = 0.5
XX=range(N)
for i in range(N):
    XX[i]=[data_s['length'][i], data_s[i_metric][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
XX = np.array(XX)

# XX[ length, i , volume, max_i]
loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))

cso_m = CSO_it[0].astype(int)

s1 = 40

qcso= cutoff
HN = tn*qcso*DR/bern_h
Hn=HN

env_y = []
Tr = np.arange(Hn+0.001,6,0.01)

for Tr1 in Tr:
    zz = -(np.exp(-Tr1/Hn)*Tr1/Hn)
#    print(zz)
    w = lambertw(zz,0)
    Tn1 = (Hn*Tr1)/(Tr1+Hn*w.real)
    env_y.append(Tn1)
#plt.plot(Tr, env_y)

tn = 3*5/3


# show non CSO storm events
fig10 = plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
ax = plt.gca()
plt.plot([HN,10],HN*np.ones(2), c='b',label = 'GN-Lower Bound')
plt.plot(HN*np.ones(2),[HN,10], c ='b')
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
#color1 =['b','y','g','r','c','k','m','g','r','lightgreen']
#Alpha = 0.5
#XX=range(N)
#for i in range(N):
#    XX[i]=[data_s['length'][i], data_s[i_metric][i],data_s['v'][i],data_s['max_i'][i],data_s['storm_id'][i]] #  duration, mean i, volume, max_i
#XX = np.array(XX)
#
## XX[ length, i , volume, max_i]
#loc_irr = list(set(CSO_re[0]).intersection(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4]))
#loc_rel = list(set(CSO_re[0]).union(CSO_re[1],CSO_re[2],CSO_re[3],CSO_re[4])- set(loc_irr)) 
#CSOs = list(set(CSO_it[0].astype(int))-set(loc_rel)- set(loc_irr))
#nonCSO=list(set(range(len(XX))) - set(CSOs) -set(loc_rel)- set(loc_irr))
#
#cso_m = CSO_it[0].astype(int)
#
#s1 = 40  # the circle size on the figure
#qcso= cutoff
#HN = tn*qcso*DR/bern_h
#a = 3/5 # adjustion fraction for the storage increase from infiltration of long duration storms
## only appied to the horozontal bound 
#env_y = []
#Hn_x=HN
#Hn_y=HN*a
#Hn=HN
#Tr = np.arange(Hn+0.001,6,0.01)
#
#for Tr1 in Tr:
#    zz = -(np.exp(-Tr1/Hn)*Tr1/Hn)
##    print(zz)
#    w = lambertw(zz,0)
#    Tn1 = (Hn*Tr1)/(Tr1+Hn*w.real)
#    env_y.append(Tn1-HN*(1-a))
##plt.plot(env_y)
## show non CSO storm events
#plt.figure(figsize=(8, 6))  # Volume vs Maximum Intensity
#ax = plt.gca()
#plt.plot([Hn_x,10],(Hn_y)*np.ones(2), c='b',label = 'GN-Adjusted Lower Bound')   # y = HN
#plt.plot(Hn_x*np.ones(2),[Hn_y,10], c ='b')                    # x = HN
#plt.scatter(XX[nonCSO,2]*DR/(bern_h), XX[nonCSO,1]*tn*DR/(bern_h),s=s1, c='grey',
#                   marker=markers1[0], alpha=0.2, edgecolor="grey", label = 'Non-CSOs')
## show CSO storms
#plt.scatter(XX[CSOs,2]*DR/(bern_h), XX[CSOs,1]*tn*DR/(bern_h),s=s1, c='k',
#                   marker=markers1[0], alpha=0.9, edgecolor="k", label = 'Persisting CSOs')
## plot Loc-irrelevent storms
#plt.scatter(XX[loc_irr,2]*DR/(bern_h), XX[loc_irr,1]*tn*DR/(bern_h), c='b',linewidths=0, marker=markers1[0], s=s1, alpha=1 ,edgecolor="b",
#                label = 'Loc-irrelevant CSOs')  
#plt.scatter(XX[loc_rel,2]*DR/(bern_h), XX[loc_rel,1]*tn*DR/(bern_h),  c='r',linewidths=0, marker=markers1[0], s=s1, alpha =Alpha,edgecolor="r",
#                label = 'Loc-relevant CSOs')
#plt.plot(Tr,env_y,label = 'GN-Ajusted Envelop',c='g')
#ax.set_yscale('log')
#ax.set_xscale('log')
#
#plt.ylabel("Tn")
#plt.xlabel("Tr")
#plt.xlim(0.01,10)
#plt.ylim(0.01,10)
#plt.rcParams.update({'font.size': 16})
#plt.legend()
#plt.legend(loc="upper left", scatterpoints=1, fontsize=12)
#
#plt.tight_layout()
#
#plt.grid()
#plt.savefig("Fig11_GN_%s_tn_%s_Tn_Tr_dr_%s_adjusted_envelope_line_storage_adjusted.png" %(city,tn,DR))


