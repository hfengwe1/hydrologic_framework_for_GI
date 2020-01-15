"""
Created on Sun May 21 20:42:17 2017
@author: Fengwei Hung

This code calculates the volumes and indices of the CSO.
The output is save in excel spreadsheet where the cso storm indices and volumes are saved by the location of GI
This will be the input of the code for making Figure 14 and 15
"""
from __future__ import division
import os
os.chdir('C:/School/2nd Paper/model')  # for my computer only
import numpy as np
from cat_runoff_w_infil import sub_runoff
from cat_runoff_w_infil import Muskingun
import pandas as pd
from def_storm_CDO import storm_id

#from Hydrograph_plot import Muskingun_plot 
cutoff = 0.20 # Q_CSO: cfs
city = 'Philadelphia'
#city = 'Seattle'
filename = city +'_1980_2013.csv'  # unit: =in/hr

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

P_crit = 0.01 # in/hr
P_lag = 3  # P_lag hours of no precipitation is consider a different storm
metric = 'in' # indicate the metric in the input file: mm or in

#%%
End_year = 2013
Start_year = 1980
run_dates =     [str(Start_year)+'-01-01', str(End_year)+'-12-31']
n_yr = End_year-Start_year
df_P1, data_s = storm_id(filename, run_dates,P_crit, P_lag,metric)   # mean_i, length, max_i, storm_id, volume(mm)
df_P = df_P1.resample('10min').bfill()
rain= df_P['P'].values.tolist()


#%% 
################################
#### CSO volume and indices ####
################################
tn = 3
cso_id =np.ones(len(data_s['storm_id']))*-1
data_s['cso_id'] = cso_id 
CSO_it = []  ## index of the storms that generate CSO
CSO_re = []  ## CSO event reduction
CSO_t=[]
Peak_t=[]
CSO_f=[]
CSO_V = []  # CSO volume
err_v = 0.001  # the monimum flow exceeding q_CSO (in) for a CSO event/ used for aoiding numberical errors
F_V = []    # cumulative infiltration

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
        for j in range(len(cso_ii)):
            data_s['cso_id'][cso_ii[j]] = 1
        
    cso_kgal=sum(np.array(cso_v))  # k-gal
    CSO_t.append(cso_kgal)
    peak = max(Q_t)
    Peak_t.append(peak)
        
#data_cso = np.zeros(len(cso_id))   
data_s['CSO_V'] = np.ones(len(data_s['storm_id']))*0
data_s['Infil'] = np.ones(len(data_s['storm_id']))*0
data_s.iloc[CSO_it[0],8] = CSO_V[0] # CSO volume without GI
data_s.iloc[CSO_it[4],9] = F_V[4]   # infiltration of the persisting CSOs


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
writer.save()
   
