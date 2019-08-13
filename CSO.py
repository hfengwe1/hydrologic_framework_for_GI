# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:19:36 2017

@author: fengwei
"""

from __future__ import division
import os
os.chdir('C:/School/2nd Paper/model')  # for my computer only
import numpy as np
from cat_runoff_w_infil import sub_runoff
from cat_runoff_w_infil import Muskingun

def cso_Mgal(tr,k,x,cutoff,inflow,rain,all_area,tot_area,imperv_perc,gi_s_p,dr,Sf,L,dt,End_T,n_sub,location,manning_n,bern_h):
    Q = inflow
    Q_all = []

    for j in range(1,(n_sub+1)):  # calculate hydrograph at the each subcatchment
        if location == j:
            gi_s = all_area*gi_s_p     # Surface area of GI (ft^2)   
            gi_d = gi_s*dr             # drainage area (ft^2)  
            w = tot_area[j-1]/L          # width (ft)
            w_gi = gi_d/L                # width
            rf_s1, f_s1, d_s1 =  sub_runoff(rain,tot_area[j-1],imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h) # [ p, imp, gi]
            rf = rf_s1[0]+ rf_s1[1]+ rf_s1[2]
            
            c0=(-k*x+0.5*dt)/(k-k*x+0.5*dt)
            c1=(k*x+0.5*dt)/(k-k*x+0.5*dt)
            c2=(k-k*x-0.5*dt)/(k-k*x+0.5*dt)
            Q = Muskingun(Q, rf,dt,End_T,c0,c1,c2,rain)
            Q_t=np.array(Q)
            Q_all.append(Q_t)
            F =np.sum(np.array(f_s1[1])*dt)/12*gi_s*7.48/1000000  # million gallon
        else:
            gi_s = 0     # Surface area of GI (ft^2)   
            gi_d = 0     # drainage area (ft^2)  
            w = tot_area[j-1]/L          # width (ft)
            w_gi = gi_d/L                # width
            rf_s1, f_s1, d_s1 = sub_runoff(rain,tot_area[j-1],imperv_perc,gi_s,gi_d,Sf,w,w_gi,dt,End_T,manning_n,bern_h) # [ p, imp, gi]
            rf = rf_s1[0]+ rf_s1[1]+ rf_s1[2]

            c0=(-k*x+0.5*dt)/(k-k*x+0.5*dt)
            c1=(k*x+0.5*dt)/(k-k*x+0.5*dt)
            c2=(k-k*x-0.5*dt)/(k-k*x+0.5*dt)
            Q = Muskingun(Q, rf,dt,End_T,c0,c1,c2,rain)  # update Q
            Q_t=np.array(Q)
            Q_all.append(Q_t)  # outlfow = Q[4].
            F =np.sum(np.array(f_s1[1])*dt)/12*gi_s*7.48/1000000  # million gallon
#### CSO Estimation 
    overflow=[ j for j in Q if j > (cutoff)]
    cso=(sum(overflow)-len(overflow)*cutoff)*dt*3600  # cubic feet
    cso_Mgal =cso*7.48/(1000000)   # M-gal
    return Q_t,Q_all,cso_Mgal,F
            