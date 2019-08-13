# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 23:10:53 2017
@author: fengwei
"""
from __future__ import division
import os
os.chdir('C:/School/2nd Paper/model')  # for my computer only
import numpy as np
from SWMM_infil import perv
from RK41 import RK4

# rain (in/hr)
# ds, depression storage (in)
# Dtheta, soil moisture deficit (unitless)
# Ks, saturated infiltration rate
# psi, suction head (in)
# eva, evaporation, (in/hr)
# n_p, manning's coefficient
# Sf, surface slope (unitless)
# w, width (ft)
# Area_p, pervious surface area (ft^2)

def sub_perv(rain,ds,Dtheta,Ks,psi,dt,eva,n_p,Sf_p,w,Area_p):
    F=0              # cumulative infiltration (in)
    f=0
    T=0              # recovery time remaining before the next event can begin 
    t=0
    d=0              # depth of ponding (in)
    q=0              # runoff flow rate (in/hr)
    Du=Dtheta        # soil moisture deficit in the upper soil recovery zone
    Dmax=0.432      # maximum soil moisture deficit 
    f_p_t=[]         # infiltration rate, f(t)  (in/hr)
    d_p_t=[]         # sotrage depth (in/hr)
    q_p_t=[]         # runoff flow rate 
    count = 0
    while count < len(rain):
        i=rain[count]
        (f,F,T,Du,da) = perv(i,d,dt,Dtheta,Du,Dmax,Ks,psi,F,T)
#        print "infiltration rate=",f# in/hr
#        print "Cumulative Infiltration =",  F # in
#        print "Soil Moisture Deficit =", Du 
#        print "t = ", t # hr
#        print " d", d  # depth of ponding
#        dd = change in water depth
        f_p_t.append(f)     # infiltration time serie (in/hr)
        # da: effective rainfall
        if i > 0:            # raining - assuming no evaporation 
            ix = i - f      # no inflow rate
            if da <= ds:
                d = max(da,0)
                q = 0
            elif d < ds and da >= ds:
#                print "Case !"
                dtp = (ds - d)/ix                             # time to ponding
                dtx = dt - dtp                               # remaing time after ponding
                tx = t +  dtp                                # time of ponding
                alpha = 1.49*w*Sf_p**0.5/Area_p/n_p    # q = alpha*d**(5/3))  (ft/s)
                dx = da - ds
                ddx = RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                # ddx: decrease in water depth
                t1=tx
                t2=tx+dtx
                while t1 <= t2:                                       # use Runge-Kutta to update dx
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)      
                    if dx <= 0:
                        dx = 0
                        break
                dx = max(dx,0)
                q = (da - (ds + dx))/dt       # runoff flow rate in/hr   
                q = max(q,0)
                d = dx + ds 
                d = max(d,0)           

            elif d > ds and da < ds:     # reception period 0 might need a second look
                dtp = (ds - da)/ix                             # time to no ponding
                dtx = dt - dtp                               # remaing time after no ponding
                tx = t +  dtp                                # time of ponding
                alpha = 1.49*w*Sf_p**0.5/Area_p/n_p    # q = alpha*d**(5/3))   (ft/s/ft^2)
                dx = d - ds                          # calculate the dx of the transition from ponding to no ponding
                ddx = RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=tx
                t2=tx+dtp
                while t1 <= t2:                                       # use Runge-Kutta to update dx
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)      
                    if dx <= 0:
                        dx = 0
                        break
                q = (d-(ds+dx))/dtp       # runoff flow rate in/hr   
                q = max(q,0)
                d = da - q*dt              # new depth = da (effective depth) - outflow  
                d = max(d,0)           
           
            else:   # d, da > ds 
                alpha = 1.49*w*Sf_p**0.5/Area_p/n_p   # q = alpha*d**(5/3))  
                dx= da -ds
                ddx= RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=t
                t2=t+dt
                h=float(10)
                while t1 < t2:
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h) 
                    if dx <= 0:
                        dx = 0
                        break
                q = (da - (ds + dx))/dt       # in/hr   
                q = max(q,0)
                d = dx + ds 
                d = max(d,0)

        else:   # i = 0  drying; no rain 
            if da <= ds:
                d = da - eva*dt
                d = max(d,0)
    #            print "case 2.3.1"    
                q = 0         
            else:   # da > ds 
                alpha = 1.49*w*Sf_p**0.5/Area_p/n_p   # q = alpha*d**(5/3))  
                dx= da - ds                 
                ddx= RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=t
                t2=t+dt
                while t1 <= t2:
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)      
                    if dx <= 0:
                        dx = 0
                        break
    #            dx += ddx(t*3600,dx,60)*12     # dx = d- ds, -> d = dx +ds  # in
                dx = max(dx,0)
                q = (da - (ds + dx))/dt       # in/hr   
                q = max(q,0)
                d = dx + ds - eva*dt
                d = max(d,0)  
        if d< 0:
            print "!! d<0 !! ", d 
        d_p_t.append(d)
        q_p_t.append(q)
        count += 1
        t += dt
#        print q,d
     
    f_p=np.array(f_p_t)*dt      # in
    q_p=np.array(q_p_t)*dt      # in
    d_p=np.array(d_p_t)         # in
    rf_p = q_p/12*Area_p     # ft^3
    return rf_p, f_p, q_p, d_p

##  ----------------------------------------------------------------------- ##
### Runoff Generated from impervious Area
##
def sub_imperv(rain,dt, ds,eva,n_imp,Sf_imp,w,Area_imp):
    d=0
    q_imp_t=[]
    d_imp_t=[]
    t = 0
    q=0              # runoff flow rate (in/hr)
    count=0
    while count < len(rain):
        i=rain[count]
        da = d + i*dt
        if i > 0:         
            if da <= ds:
                d = da
                q = 0
            elif d < ds and da >= ds:
                dtp = (ds - d)/i                            # time to ponding
                dtx = dt - dtp                         # remaing time after ponding
                tx = t +  dtp                                # time of ponding
                alpha = 1.49*w*Sf_imp**0.5/Area_imp/n_imp   # q = alpha*d**(5/3))  
                dx = da - ds
                ddx = RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=tx
                t2=tx+dtx
                while t1 <= t2:                                       # use Runge-Kutta to update dx
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)
                    if dx <= 0:
                        dx = 0
                        break
                dx = max(dx,0)
                q = (da - (ds + dx))/dt       # runoff flow rate in/hr   
                q = max(q,0)
                d = dx + ds 
                d = max(d,0)             
            else:   # da > ds 
                alpha = 1.49*w*Sf_imp**0.5/Area_imp/n_imp   # q = alpha*d**(5/3))  
                dx= da - ds
                ddx= RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=t
                t2=t+dt
                while t1 <= t2:
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)
                    if dx <= 0:
                        dx = 0
                        break                    
                dx = max(dx,0)
                q = (da - (ds + dx))/dt       # in/hr   
                q = max(q,0)
                d = dx + ds 
                d = max(d,0) 
        else:   # i = 0  drying
            if d <= ds:
                d += - eva*dt
                d = max(d,0)
                q = 0         
            else:   # da > ds 
                alpha = 1.49*w*Sf_imp**0.5/Area_imp/n_imp   # q = alpha*d**(5/3))  
                dx=d - ds 
                ddx= RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=t
                t2=t+dt
                while t1 <= t2:
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)          
                    if dx <= 0:
                        dx = 0
                        break
                dx = max(dx,0)
                q = (da - (ds + dx))/dt       # in/hr   
                q = max(q,0)
                d = dx + ds - eva*dt
                d = max(d,0)       
### Check Mass Balance ####
#        if i> 0 :
#            mass_err = d -(da - q*dt)
#        elif i==0 and d > 0:
#            mass_err = d -(da - q*dt-eva*dt) 
#        else:
#            mass_err = d -(da - q*dt) 
#            
#        print mass_err
        if np.isnan(d):
            print'd = nan'
            break
        
        d_imp_t.append(d)
        q_imp_t.append(q)
        count += 1
        t += dt
     
    q_imp=np.array(q_imp_t)  # in/hr
    d_imp=np.array(d_imp_t)  # in/hr
#    Rain = np.array(rain)*dt
#    sum(Rain)
#    sum(data_q_imp)*dt
    rf_imp = (q_imp/12)*dt*Area_imp     # ft^3
    
#### Check Mass Balance ####
#mass_err = sum(np.array(rain)*dt/12*Area_imp) - sum(rf_imp) - eva/12*dt*(gi_d-gi_s)*8 
# assumed that 8*dt hours are drying (eva)

    return rf_imp,q_imp,d_imp
##########################
#### GI: rain garden
##########################

# bern_h: in
# s_dep: soil depth (in), s_Dtheta: soil porosity, s_Ks: saturated conductivity
# s_psi: suction head (in)     
def sub_gi(rain,Dtheta,Ks,psi,dt,gi_d,gi_s,bern_h,Area_gi,w_gi,n_gi,Sf_gi):
    F=0              # cumulative infiltration (in)
    f=0
    T=0 
    t=0
    d=0              # initial depth of ponding (in)
#    tot_infil=0      # total infiltration (in)
    Du=Dtheta        # soil moisture deficit in the upper soil recovery zone
#    Dmax=0.432      # maximum soil moisture deficit   
    Dmax=0.001      # maximum soil moisture deficit   
    
    f_gi_t=[]           # infiltration rate, f(t) (in/hr)
    d_gi_t=[]
    of_gi_t=[]          # runoff geneation rate (in/hr)
#    End_T = (len(rain))*dt
    count = 0
  
#### Runoff Generated from Impervious Area ####
    ds=0.00             # depression storage (in)
    eva = 0.05         # evaporation rate (in/hr)
    Area_imp=gi_d-gi_s
    n_imp = 0.015    # Manning Coefficient
    Sf_imp = 0.03   # % 
    if Area_imp>0:    
        rf_imp,q_imp,d_imp = sub_imperv(rain,dt, ds,eva,n_imp,Sf_imp,w_gi,Area_imp)
    else:
        rf_imp = np.zeros(len(rain))
#        q_imp = rf_imp
#        d_imp = rf_imp
    # rf_imp: ft^3
    # q_imp: in/hr
    # d_imp: in
    ds = bern_h # bern_h = bern hight (in)

    while count < len(rain):
        i=rain[count]+rf_imp[count]*12/gi_s/dt    # rainfall intensity (in/hr)
        (f,F,T,Du,da) = perv(i,d,dt,Dtheta,Du,Dmax,Ks,psi,F,T)  # calculate infiltration rate
#        if d >= bern_h:
#            of=(d-bern_h)/12*Area_gi  # ft^3
#            d=bern_h
#        else:
#            of=0
        if i > 0:            # raining - assuming no evaporation 
            ix = i - f      # effective i 
            if da <= ds:
                d = max(da,0)
                q = 0
            elif d < ds and da >= ds:
#                print "Case !"
                dtp = (ds - d)/ix                             # time to ponding
                dtx = dt - dtp                               # remaing time after ponding
                tx = t +  dtp                                # time of ponding
                alpha = 1.49*w_gi*Sf_gi**0.5/Area_gi/n_gi    # q = alpha*d**(5/3))  (ft/s)
                dx = da - ds
                ddx = RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                # ddx: decrease in water depth
                t1=tx
                t2=tx+dtx
                while t1 <= t2:                                       # use Runge-Kutta to update dx
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)      
                    if dx <= 0:
                        dx = 0
                        break
                dx = max(dx,0)
                q = (da - (ds + dx))/dt       # runoff flow rate in/hr   
                q = max(q,0)
                d = dx + ds 
                d = max(d,0)           

            elif d > ds and da < ds:     # reception period 0 might need a second look
                dtp = (ds - da)/ix                             # time to no ponding
                dtx = dt - dtp                               # remaing time after no ponding
                tx = t +  dtp                                # time of ponding
                alpha = 1.49*w_gi*Sf_gi**0.5/Area_gi/n_gi    # q = alpha*d**(5/3))   (ft/s)
                dx = d - ds                          # calculate the dx of the transition from ponding to no ponding
                ddx = RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=tx
                t2=tx+dtp
                while t1 <= t2:                                       # use Runge-Kutta to update dx
                    h=float(10)
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h)      
                    if dx <= 0:
                        dx = 0
                        break
                q = (d-(ds+dx))/dtp       # runoff flow rate in/hr   
                q = max(q,0)
                d = da - q*dt              # new depth = da (effective depth) - outflow  
                d = max(d,0)           
           
            else:   # d, da > ds 
                alpha = 1.49*w_gi*Sf_gi**0.5/Area_gi/n_gi   # q = alpha*d**(5/3))  
                dx= da -ds
                ddx= RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=t
                t2=t+dt
                h=float(10)
                while t1 < t2:
                    t1, dx = t1+h/3600, dx-ddx(t1*3600,dx,h) 
                    if dx <= 0:
                        dx = 0
                        break
                q = (da - (ds + dx))/dt       # in/hr   
                q = max(q,0)
                d = dx + ds 
                d = max(d,0)

        else:   # i = 0  drying; no rain 
            if da <= ds:
                d = da - eva*dt
                d = max(d,0)
    #            print "case 2.3.1"    
                q = 0         
            else:   # d > ds 
                alpha = 1.49*w_gi*Sf_gi**0.5/Area_gi/n_gi     # q = alpha*d**(5/3))  
                dx= da - ds                 
                ddx= RK4(lambda t, dy:alpha*(dy/12)**(5/3))  # d(d)/dt = i - f- d - q # t = t+dt
                t1=t
                t2=t+dt
                while t1 <= t2:
                    h=float(10)
                    t1, dx = t1+h/3600, dx+ddx(t1*3600,dx,h)      
                    if dx <= 0:
                        dx = 0
                        break
                dx = max(dx,0)
                q = (da - (ds + dx))/dt       # in/hr   
                q = max(q,0)
                d = dx + ds - eva*dt
                d = max(d,0)  
         
    
#        print " ----------------------Separate line------------------------"
#        print "t=", count
#        print "infiltration rate=",f
#        print "Cumulative Infiltration =", F
#    #    print "Water-in ", sum(np.array(rain[0:(count+1)])*dt)*gi_d/gi_s, " = ", "Water-out " , tot_infil+d+sum(np.array(rf_t[0:(count+1)]))
#        print "Soil Moisture Deficit =", Du 
#        print "Storage = ", d, "in"
#        print "Overflow = ", q, "in/hr"

#######  Check Mass Balance   ########
#        if i > 0: 
#            mass_err= i*dt + (d_gi_t[count] -d) - (q + f)*dt
#        else: 
#            mass_err= i*dt + (d_gi_t[count] -d) - (q + f)*dt - eva*dt
#        print mass_err
            
        d_gi_t.append(d)
        f_gi_t.append(f)    # infiltration rate (in/hr)
        of_gi_t.append(q)  # overflow  (in/hr)
        count += 1
        t += dt
#        print (i,f,F,d,of)
#        print " t = ? ", t  
    f_gi=np.array(f_gi_t)              # in/hr
    d_gi=np.array(d_gi_t)              # in
    of_gi=np.array(of_gi_t)/12*dt*gi_s # ft^3
 
    return of_gi, f_gi, d_gi, 
#### Check Mass Balance  ####   
#    d_end= d_gi[len(d_gi)-1]
#    storage1= d_end*gi_s/12
#    infiltration1 = sum( f_gi/12*dt*gi_s)
#    evaporation1=rain.count(0)*eva/12*dt*gi_s
#    of_mass= sum(of_gi)
##    dep_s= ds*(gi_d-gi_s)/12
#    mass_out= storage1+infiltration1 + evaporation1 + of_mass 
#    mass_in = sum(np.array(rain)/12*dt*gi_d)
    
#


