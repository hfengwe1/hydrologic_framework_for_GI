# -*- coding: utf-8 -*-
"""
Created on Mon Apr 03 22:55:08 2017

@author: fengwei
"""
from __future__ import division
import math as ma
import scipy.optimize
#import numpy as np
#from Routing import Manning
## SWMM Infiltration Computation Scheme
# i: rainfall intensity(in/hr), d: depth of ponding (in)
# Dtheta: soil moisture deficit, Du: # soil moisture deficit in the upper soil recovery zone
# Dmax: maximum soil moisture deficit.
# Ks: Saturated hydraulic conductivity (in/hr), psi: suction head (in), dt: time step (hr) 
# Assume that ponding depth is negligible
def perv(i,d,dt,Dtheta,Du,Dmax,Ks,psi,F,T):
    Lu=4*Ks**0.5      # depth of upper soil recovery zone
#    Tr=4.5/(Ks**0.5)  # minimum recovery time remaining (hr)
    Tr = 24
    ia=i+d
    f=0               # in/hr
    if ia == 0:      # when drying (rain stopped)
        f=0
#        theta_dt=kr*Dmax*dt
        dF= -(Dmax-Du)*dt*Lu               # drying take 1 hr
        Du = min(Du+(Dmax-Du)/4, Dmax)     # update soil moisture deficit, but not larger Dmax

        F += dF   # update cumulative infil. (where does the water go?)
        T += -dt
        if T <= 0:             # reset Dtheta and F
            Du=Dmax
            F=0
            T=0
    elif ia <= Ks:   # unsaturated condition - low intensity
        f= ia
        dF=ia*dt
        F += dF
        Du = max(Du-dF/Lu,0)
    elif ia> Ks:     # unsaturated condition - mid-high intensity
        T=Tr
        Fs=Ks*psi*Dtheta/(ia-Ks)
        if F >= Fs:  # saturated condition
            F1=F
            C=Ks*dt+F1-psi*Dtheta*ma.log(F1+psi*Dtheta)
            F2=scipy.optimize.fsolve(lambda x: C+psi*Dtheta*ma.log(x+psi*Dtheta)-x ,0.5)[0] 
            dF=F2-F
            if dF > ia*dt:
                dF=ia*dt
                F += dF
                Du = max(Du-dF/Lu,0)
                f = dF/dt
            else:
                F += dF
                Du = max(Du-dF/Lu,0)
                f = dF/dt
        elif F+ia*dt < Fs:
            f =ia
            dF=ia*dt
            F += dF
            Du = max(Du-dF/Lu,0)
        else:   # F+ia*dt >Fs and F < Fs
            F1=Fs
            dt1=dt-(Fs-F)/ia        # time of saturation
# calculate the portion of the original time step over which saturated conditions exist
            C=Ks*dt1+F1-psi*Dtheta*ma.log(F1+psi*Dtheta)
            F2=scipy.optimize.fsolve(lambda x: C+psi*Dtheta*ma.log(x+psi*Dtheta)-x ,0.5)[0] 
            dF=F2-F
            F += dF
            Du = max(Du-dF/Lu,0)
            f = dF/dt
    dd = (i - f)*dt
#    print "check mass ballance (i-f)*dt-dd = "
    d += dd  # ia = i +d/dt)
#    print f, F

    return f,F,T,Du,d

#   f: infiltration rate, in/hr
#   F: cumuative infiltration, in
#   T: timer to calculate the soil moisture
        
        
        
        
        
        


