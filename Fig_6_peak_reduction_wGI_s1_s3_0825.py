
"""
Created on Sun May 21 20:42:17 2017

@author: Fengwei Hung
"""
from __future__ import division
import os
os.chdir('C:/School/2nd Paper/model')  # for my computer only
import numpy as np
import matplotlib.pyplot as plt
from cat_runoff_no_infil import sub_runoff
from cat_runoff_no_infil import Muskingun
#from cat_runoff_w_infil import sub_runoff
#from cat_runoff_w_infil import Muskingun
from def_storm_CDO import storm_id

#from Hydrograph_plot import Muskingun_plot 
#
#tot_rain = 4 
dt =1/6

#Tr = 1
#Tn = 0.2
#tr = 3
#tgi = tr/Tr
#tn =Tn*tgi
##tgi = 2
##tr = Tr*tgi # rainfall duration,  hr
##rain_i = V/tr # rainfall intensity, in/hr
#rain_i = 0.5 # in/hr
#V = tr*rain_i # in of rainfall
#rain= list(np.ones(int(tr/dt))*rain_i)
#rain= rain +list(np.zeros(int(8/dt)))
#bern_h= 1.5
#DR=bern_h/(rain_i*tgi)
#gi_s_p = 0.05


#%% Change tr: duration and tn: travel time
rain_i = 0.5 # in/hr

DR=6
bern_h= 6
gi_s_p = 0.05
Tr = 5
Tn = 3
tgi = bern_h/(rain_i*DR)
tn =Tn*tgi
tr = Tr*tgi

#tgi = 2
#tr = Tr*tgi # rainfall duration,  hr
#rain_i = V/tr # rainfall intensity, in/hr

V = tr*rain_i # in of rainfall
rain= list(np.ones(int(tr/dt))*rain_i)
rain= rain +list(np.zeros(int(15/dt)))
#
#

#%%

End_T= (len(rain))*dt
tot_area = 43560*2.471  # surface area (ft^2) / 1 ac = 43560 ft^2
nsub = 3
sub_area = tot_area/nsub
imperv_perc = 1.0   # percentage of impervious area

Sf=[0.03,0.03,0.03]   # slope % [s_p, s_imp, s_gi]
manning_n=[0.01, 0.01, 0.01]
L=100             # length of the watershed

location = 3
inflow = np.zeros(len(rain))
k=Tn*tgi/nsub## hr parameter for Muskingun  

x=0.1 
# subcatchment discharge without GI
gi_s=0    # no GI # Surface area of GI (ft^2)  
gi_d = gi_s*DR     # drainage area (ft^2)  
w=sub_area/L      # width of the watershed
w_gi =gi_d/L
rf_s1, f_s1, d_s1 = sub_runoff(rain,sub_area,imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h) # [ p, imp, gi]
rf = rf_s1[0]+rf_s1[1]+rf_s1[2]
rf11=rf.tolist()
rf11.insert(0,0)
rf1=np.array(rf11)    
if len(rf1)>len(rain):
    rf1 = rf1[0:len(rain)]
    
# subcatchment discharge with GI
gi_s = tot_area*gi_s_p    # Surface area of GI (ft^2)
gi_d = gi_s*DR     # drainage area (ft^2)  
w=sub_area/L      # width of the watershed
w_gi =gi_d/L
rf_s3, f_s3, d_s3 = sub_runoff(rain,sub_area,imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h)        
rf = rf_s3[0]+rf_s3[1]+rf_s3[2]
rf11=rf.tolist()
rf11.insert(0,0)
rf2=np.array(rf11)
if len(rf2)>len(rain):
    rf2 = rf2[0:len(rain)]
#### Muskingun Routing

c0=(-k*x+0.5*dt)/(k-k*x+0.5*dt)
c1=(k*x+0.5*dt)/(k-k*x+0.5*dt)
c2=(k-k*x-0.5*dt)/(k-k*x+0.5*dt)
## calculate discharge without GI  
inflow1 =np.zeros(len(rf1))
   

for i in np.arange(1,(1+nsub),1):
    Q1 = Muskingun(inflow1, rf1,dt,End_T,c0,c1,c2,rain)  # no GI
    inflow1 = Q1

# Gi at S1 and S3
Q2 = Muskingun(inflow, rf2,dt,End_T,c0,c1,c2,rain)   # Gi at S1
inflow2 = Q2
Q2 = Muskingun(inflow2, rf1,dt,End_T,c0,c1,c2,rain)   # Gi at S1
inflow2 = Q2
Q2 = Muskingun(inflow2, rf1,dt,End_T,c0,c1,c2,rain)   # Gi at S1


Q3 = Muskingun(inflow, rf1,dt,End_T,c0,c1,c2,rain)   # Gi at S3, CFS
inflow3 = Q3
Q3 = Muskingun(inflow3, rf1,dt,End_T,c0,c1,c2,rain)   # Gi at S3, CFS
inflow3 = Q3
Q3 = Muskingun(inflow3, rf2,dt,End_T,c0,c1,c2,rain)   # Gi at S3 , CFS

print("Tr = %s" %Tr)    
print('Tn = %s'%Tn)
print("DR = %s" %DR)
print("tgi = %s" %tgi)
print('tn = %s'%tn)
print("tr = %s" %tr)
#%% Plot

plt.rcParams.update({'font.size': 20})

tr_1 = np.arange(0,10.5,0.5)

#qp1 = rain_i*(1- np.exp(-tr_1/1))
#qp2 = rain_i*(1- np.exp(-tr_1/tn))  # tn* = tn/2
#qp3 = rain_i*(1- np.exp(-tr_1/1.4))
#qp4 = rain_i*(1- np.exp(-tr_1/4))
#qp5 = rain_i*(1- np.exp(-tr_1/5))
#qp6 = rain_i*(1- np.exp(-tr_1/1.4))
#### SI unit:
#Qm1 = np.array(Q1)*43560/12*0.0283168466/1e4*3600   #m^3/h/ha
#Qm2 = np.array(Q2)*43560/12*0.0283168466/1e4*3600   #m^3/h/ha
#Qm3 = np.array(Q3)*43560/12*0.0283168466/1e4*3600   #m^3/h/ha
#
Qm1 = np.array(Q1)*0.0283168466*3600*1e3*1e-4   #mm/h/ha
Qm2 = np.array(Q2)*0.0283168466*3600*1e3*1e-4   #mm/h/ha
Qm3 = np.array(Q3)*0.0283168466*3600*1e3*1e-4   #mm/h/ha

fig4= plt.figure(4, figsize=(8, 5))

#x1= np.arange(0,(End_T+dt),dt)
x1= np.arange(0,(End_T),dt)

plt.plot(x1,Qm1,color="k", linewidth=1.5, linestyle="-.", label='No GI')
plt.plot(x1,Qm2,color="g", linewidth=1.5, linestyle="-", label='GI at Upper')
plt.plot(x1,Qm3,color="b", linewidth=1.5, linestyle="-", label='GI at Lower')
#plt.plot(tr_1, qp1, linewidth=1.5, linestyle="-", label='tn=1')
#
#plt.plot(tr_1, qp2, linewidth=1.5, linestyle="-", label='qp=i*(1-exp(-tr/tn)')
#plt.plot(tr_1, qp3, linewidth=1.5, linestyle="-",Label ='tn=1.4')
#plt.plot(tr_1, qp4, linewidth=1.5, linestyle="-",Label ='tn=4')
#plt.plot(tr_1, qp5, linewidth=1.5, linestyle="-",Label ='tn=5')
#plt.plot(tr_1, qp6, linewidth=1.5, linestyle="-", label='tn=1.4')


plt.ylim(0,20)
plt.xlim(0,25)
plt.rcParams['font.size'] = 20
plt.rcParams['legend.fontsize'] = 'medium'
plt.rcParams['figure.titlesize'] = 'medium'
plt.xlabel("Time (h)")
plt.ylabel("Discharge (mm/h)")
#plt.legend(loc='upper right')
plt.tight_layout()

plt.savefig("Fig_6_Tr_5_Tn_3_dr_6_tgi_2" )
#



