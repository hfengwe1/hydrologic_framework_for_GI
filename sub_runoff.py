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
import subcatchment 

def sub_runoff(rain,tot_area,imperv_perc,gi_s,gi_d,Sf,w,dt,End_T):
    #### Runoff Generated from Pervious Area
    ds = 0.05           # depression storage (in)
    Dtheta=0.2     # soil moisture deficit 
    Ks=0.05          # saturated hydraulic conductivity (in/hr)
    psi=10.75        # suction head (in)
    eva = 0.05       # evaporation rate (in/hr)
    n_p=0.03    # Manning coefficient for short grass
    area= [tot_area*imperv_perc-gi_d,tot_area*( 1-imperv_perc),gi_d]  # area of imperv. and perv. surface
    Area_p=area[1]
    x= np.arange(0,End_T,dt)
    
    rf_p, f_p, q_p, d_p = sub_perv(rain,ds,Dtheta,Ks,psi,dt,eva,n_p,Sf,w,Area_p)
    # rf_p: total volume of runoff (ft3), 
    # f_p: infiltration rate (in.hr)
    # d_p: water depth at surface (in)
    # q_p: runoff flow rate (in/hr)
#    rf_p=np.zeros(len(x))
#    rf_p[0:len(rf)]=rf    
    #### Runoff Generated from Impervious Area ####
    ds=0.05             # depression storage (in)
    eva = 0.01         # evaporation rate (in/hr)
    Area_imp=area[0]
    n_imp = 0.015    # Manning Coefficient
    Sf_imp = 0.003   # %
    
    rf_imp,q_imp,d_imp = sub_imperv(rain,dt, ds,eva,n_imp,Sf_imp,w,Area_imp)
#    rf_imp=np.zeros(len(x))
#    rf_imp[0:len(rf_i)]=rf_i 
    # rf_imp: total volume of overflow (ft^3), 
    # q_imp_t: runoff flow rate (in/hr)
    # d_imp_t: water depth at surface (in)
    
    #### Overflow from GI #########################
    if gi_d >= imperv_perc*tot_area:
        gi_d= imperv_perc*tot_area
        print "GI treats all the impervious area", gi_d, "ft^2"
    bern_h=6.0         # bern hight (in)
    Dtheta=0.332     # soil moisture deficit 
    Ks=0.05          # saturated hydraulic conductivity (in/hr)
    psi=10.75        # suction head (in)
    Area_gi=area[2]
    w_gi =(gi_d-gi_s)/50  # assume the furtherest point to the edge of a GI unit is 50 ft
#    w_gi =Area_gi**0.5
    n_gi=0.03    # Manning coefficient for short grass
    Sf_gi=0.001
    if gi_s>0:
         of_gi, f_gi, d_gi = sub_gi(rain,Dtheta,Ks,psi,dt,gi_d,gi_s,bern_h,Area_gi,w_gi,n_gi,Sf_gi)
    else:
        f_gi=np.zeros(len(x))
        d_gi=np.zeros(len(x))
        of_gi=np.zeros(len(x))    
    # f_gi_t: infiltration rate (in.hr)
    # d_gi_t: water depth at GI surface (in)
    # of_gi_t: overflow rate (in/hr)    
    rf_gi= np.array(of_gi)*dt      #ft^3      
    sub_rf = [rf_p, rf_imp,rf_gi]  #ft3
    sub_f=[f_p, f_gi]
    sub_d=[d_p, d_imp, d_gi]    
    return sub_rf, sub_f, sub_d
    
def kinematicWave(rf_out_t,rain,dt,n_s,Bs,L,s0,End_T):
#### Routing - Kinematic wave
    ####  outflow (ft^3/s)
    #### Routing  from sub 1 to sub 2
    t = 0

    rf_p= rf_out_t[0]
    rf_imp= rf_out_t[1]
    rf_gi= rf_out_t[2]
    Q1= np.zeros(int(End_T/dt))
    Q1[0:len(rain)] = (rf_p + rf_imp + rf_gi)/(dt*3600)  # CFS
    np.insert(Q1,0,[0])   # insert Q = 0 at t = 0
    y_t = []
    ck_t = []
    dt2_t = []
    t2_t = []
    t2=0
    count = 0
    while t < End_T:
        if Q1[count] >0:
            y = (Q1[count]*n_s/(1.49*Bs*s0**0.5))**(0.6)   # y = (d-ds)
            q =1.49*(s0)**0.5/n_s*(y)**(2/3)
            ck = 5./3. *q  # Ck = 5/3 *q ft/sec
            alpha=1.49/n_s**0.5
            dt2 = (L/(alpha*y**(2/3)))**(3/5)        # min
            t2 = t*60 + dt2 # min
        else:
            y = 0
            ck = 0
            dt2 = 0
            t2 = t2
        y_t.append(y)
        ck_t.append(ck)
        dt2_t.append(dt2)
        t2_t.append(t2)
        t += dt
        count += 1   
    return y_t,ck_t,dt2_t,t2_t,Q1
 

def Muskingun(Q1, rf_out_t,dt,End_T,c0,c1,c2,rain):
    t=dt
    rf_p= rf_out_t[0]
    rf_imp= rf_out_t[1]
    rf_gi= rf_out_t[2]
    x= np.arange(0,End_T,dt)
    Inflow=np.zeros(len(x))
    Inflow[0:len(rf_p)] = Q1+(rf_p + rf_imp + rf_gi)/(dt*3600)  # CFS 
    Outflow = [0]
    count = 0
    while t < End_T:
        O2=c0*Inflow[count+1]+c1*Inflow[count]+c2*Outflow[count]
        Outflow.append(O2)
        t += dt
        count += 1
    return Outflow

       
    
    
    
    
