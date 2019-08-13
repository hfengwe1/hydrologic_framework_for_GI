# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 21:42:36 2017

@author: fengwei
"""
from __future__ import division
import os
os.chdir('C:/School/2nd Paper/model')  # for my computer only
import numpy as np
#import matplotlib.pyplot as plt
#from Routing import Manning
from subcatchment import sub_perv
from subcatchment import sub_imperv
from subcatchment import sub_gi

## sub_runoff- calculate runoff generates from pervious, impervious area and overflows from GI in a catchment 
def sub_runoff(rain,s_area,imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h):
    #### Runoff Generated from Pervious Area
    ds = 0.05        # depression storage (in)
    Dtheta=0.2       # soil moisture deficit 
    Ks = 0.05          # saturated hydraulic conductivity (in/hr)
    psi = 10.75        # suction head (in)
    eva = 0.05       # evaporation rate (in/hr)
    
    n_p=manning_n[0]    # Manning coefficient for short grass
    n_imp = manning_n[1]     # Manning Coefficient
    n_gi=manning_n[2]    # Manning coefficient for short grass
    
    Sf_p = Sf[0]
    Sf_imp = Sf[1]   # %  
    Sf_gi= Sf[2]
    if gi_d >= imperv_perc*s_area:
        gi_d= imperv_perc*s_area
        print "GI treats all the impervious area", gi_d, "ft^2"
    sub_area= [s_area*imperv_perc-gi_d,s_area*( 1-imperv_perc),gi_d]  # area of imperv. and perv. surface
    Area_p=sub_area[1]
    x= np.arange(0,End_T,dt)
    if Area_p > 1000:   # for area >1000 sqft
        rf_p, f_p, q_p, d_p = sub_perv(rain,ds,Dtheta,Ks,psi,dt,eva,n_p,Sf_p,w,Area_p)
    elif Area_p>10 and Area_p <=1000:  # area in between 10 and 1000 sqft 
        w_p = Area_p**0.5   
        rf_p, f_p, q_p, d_p = sub_perv(rain,ds,Dtheta,Ks,psi,dt,eva,n_p,Sf_p,w_p,Area_p)
    else:
        rf_p = np.zeros(len(rain))
        f_p = rf_p
        d_p =rf_p
# rf_p: total volume of runoff (ft3), 
# f_p: infiltration rate (in.hr)
# d_p: water depth at surface (in)
# q_p: runoff flow rate (in/hr)
#### Runoff Generated from Impervious Area ####
    ds=0.05             # depression storage (in)
    eva = 0.01          # evaporation rate (in/hr)
    Area_imp=sub_area[0]

    if Area_imp > 1000:  # if Area_imp < 100 ft^2, it will result in numberical errors       
        rf_imp,q_imp,d_imp = sub_imperv(rain,dt, ds,eva,n_imp,Sf_imp,w,Area_imp)
    elif Area_imp>10 and Area_imp <=1000: 
        w_imp = Area_imp**0.5
        rf_imp,q_imp,d_imp = sub_imperv(rain,dt, ds,eva,n_imp,Sf_imp,w_imp,Area_imp)
    else:
        rf_imp = np.zeros(len(rain))
        d_imp = np.zeros(len(rain))  
# rf_imp: total volume of overflow (ft^3), 
# q_imp_t: runoff flow rate (in/hr)
# d_imp_t: water depth at surface (in)

#### Overflow from GI #########################
    Dtheta=0.332     # soil moisture deficit 
    Ks=0.05          # saturated hydraulic conductivity (in/hr)
    psi=10.75        # suction head (in)
    Area_gi=sub_area[2]


    if gi_s>1000:
        of_gi, f_gi, d_gi = sub_gi(rain,Dtheta,Ks,psi,dt,gi_d,gi_s,bern_h,Area_gi,w_gi,n_gi,Sf_gi)
    elif gi_s>10 and gi_s <=1000: 
        w_gi2 = Area_imp**0.5
        of_gi, f_gi, d_gi = sub_gi(rain,Dtheta,Ks,psi,dt,gi_d,gi_s,bern_h,Area_gi,w_gi2,n_gi,Sf_gi)
    else:
        f_gi=np.zeros(len(x))
        d_gi=np.zeros(len(x))
        of_gi=np.zeros(len(x))    
    # f_gi_t: infiltration rate (in.hr)
    # d_gi_t: water depth at GI surface (in)
    # of_gi_t: overflow rate (in/hr)    
    rf_gi= np.array(of_gi)      #ft^3      
    sub_rf = [rf_p, rf_imp,rf_gi]  #ft3
    sub_f=[f_p, f_gi]
    sub_d=[d_p, d_imp, d_gi]    
    return sub_rf, sub_f, sub_d
    
def Muskingun(inflow, rf,dt,End_T,c0,c1,c2,rain):
    t=dt

    x1= np.arange(0,End_T,dt)
    Inflow=np.zeros(len(x1))
    Inflow[0:len(rf)] = inflow+(rf)/(dt*3600)  # CFS 
    Outflow = [(1-c1)*Inflow[0]]
    count = 0
    while count < (len(rf)-1):
        O2=c0*Inflow[count+1]+c1*Inflow[count]+c2*Outflow[count]
        Outflow.append(O2)
#        print Inflow[count], O2
        t += dt
        count += 1
    return Outflow

#def kinematicWave(rf_out_t,rain,dt,n_s,Bs,L,s0,End_T):
##### Routing - Kinematic wave
#    ####  outflow (ft^3/s)
#    #### Routing  from sub 1 to sub 2
#    t = 0
#
#    rf_p= rf_out_t[0]
#    rf_imp= rf_out_t[1]
#    rf_gi= rf_out_t[2]
#    Q1= np.zeros(int(End_T/dt))
#    Q1[0:len(rain)] = (rf_p + rf_imp + rf_gi)/(dt*3600)  # CFS
#    np.insert(Q1,0,[0])   # insert Q = 0 at t = 0
#    y_t = []
#    ck_t = []
#    dt2_t = []
#    t2_t = []
#    t2=0
#    count = 0
#    while t < End_T:
#        if Q1[count] >0:
#            y = (Q1[count]*n_s/(1.49*Bs*s0**0.5))**(0.6)   # y = (d-ds)
#            q =1.49*(s0)**0.5/n_s*(y)**(2/3)
#            ck = 5./3. *q  # Ck = 5/3 *q ft/sec
#            alpha=1.49/n_s**0.5
#            dt2 = (L/(alpha*y**(2/3)))**(3/5)        # min
#            t2 = t*60 + dt2 # min
#        else:
#            y = 0
#            ck = 0
#            dt2 = 0
#            t2 = t2
#        y_t.append(y)
#        ck_t.append(ck)
#        dt2_t.append(dt2)
#        t2_t.append(t2)
#        t += dt
#        count += 1   
#    return y_t,ck_t,dt2_t,t2_t,Q1
#        
    
    
    
    

