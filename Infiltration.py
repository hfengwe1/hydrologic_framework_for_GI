# -*- coding: utf-8 -*-
"""
Infiltration methods: Horton and Green-Ampt
Created on Mon Mar 27 10:48:04 2017
@author: fengwei
"""

import math as ma
import scipy.optimize

def horton(t,dt,fc,f0,k):
        ft=(fc+(f0-fc)*ma.exp(-k*t*dt))
        Ft=(fc*t*dt+(f0-fc)/k*(1-ma.exp(-k*t*dt)))
        return ft,Ft

def GreenAMPT(t,h_suc,Dtheta,Ks):
    Ft=scipy.optimize.fsolve(lambda x: Ks*t+h_suc*Dtheta*ma.log((x+h_suc*Dtheta)/(h_suc*Dtheta))-x ,0.5)[0] 
    ft=K*(h_suc*Dtheta/Ft +1)
    return (ft,Ft)   



#(f1,F1)=horton(fc,f0,k,T)
#(f2,F2)=GreenAMPT(t,h_suc,e_theta,porosity,K,Se)
