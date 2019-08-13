# -*- coding: utf-8 -*-

"""
Created on Sun May 21 20:42:17 2017

@author: Fengwei Hung
"""
from __future__ import division
import os
os.chdir('C:/School/2nd Paper/model')  # for my computer only
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from cat_runoff_w_infil import sub_runoff
from CSO import cso_Mgal 
np.seterr(divide='ignore', invalid='ignore')

## SI Unit 
dt =1/12
tr=2 # rainfall duration,  hr
rain_i = 0.5 # rainfall intensity, in/hr
rain= list(np.ones(int(tr/dt))*rain_i)
rain= rain + list(np.zeros(int(20/dt)))
End_T= len(rain)*dt
#all_area = 2.471*43560  # surface area (ft^2) / 1 ac = 43560 ft^2
all_area = 2.471*43560  # surface area (ft^2) / 1 ac = 43560 ft^2

n_sub=5      # number of subwatershed
sub_p= np.array([0.2,0.2,0.2,0.2,0.2])
tot_area = all_area*sub_p
imperv_perc = 1.0   # percentage of impervious area
gi_s_p= 0.025             # gi surface area to total area (percentage)
gi_s = all_area*gi_s_p     # Surface area of GI (ft^2)   (1%)

Sf=[0.03,0.03,0.03]   # slope % [s_p, s_imp, s_gi]
manning_n=[0.01, 0.01, 0.01]
DR= 6   # drainage area ratio
gi_d = gi_s*DR     # drainage area (ft^2)  
L=100            # length of the watershed
w=tot_area/L      # width of the watershed
w_gi =gi_d/L
location= 0 # lacation of GI
bern_h=6 # GI's storage depth
## Plot w/o GI
gi_s = 0   # Surface area of GI (ft^2)   
gi_d = 0     # drainage area (ft^2)  
rf_s, f_s, d_s = sub_runoff(rain,tot_area[0],imperv_perc,gi_s,gi_d,Sf,w[0],w_gi,dt,End_T,manning_n,bern_h) # [ p, imp, gi] 
rf = (rf_s[0]+rf_s[1]+rf_s[2])/(dt*3600)     # cfs   

####  Open Channel (sewer)
w_c = 100 # ft
L_c = 3000 # ft
g = 32.174 # ft^2/s
h = rf/w_c # ft
celerity = (g*h)**0.5 
ce_avg = np.average(celerity[0:int((tr/dt))*2])  # ft /s
#### Muskingun Routing

#k = L_c/ce_avg/3600  # hr
k = 0.08 # hr
x = 0.131 
#if k> 2*k*x and k<2*k*(1-x):
#    print " k value is ok! (k=%s)" %k
#else:
#    print " k value is not appropreiate: k= %s" %k
##k=0.4 ## hr parameter for Muskingun  
#
#c0=(-k*x+0.5*dt)/(k-k*x+0.5*dt)
#c1=(k*x+0.5*dt)/(k-k*x+0.5*dt)
#c2=(k-k*x-0.5*dt)/(k-k*x+0.5*dt)
#inflow=np.zeros(len(rain))

#### Timescales
#t_gi = 6/(rain_i*DR)
#t_n = k*n_sub 
#T_r=tr/t_gi
#T_n = t_n/t_gi

#### Create Tr, Tng evauation points
x1= np.arange(0.5,(dt+End_T),dt)
ix = 21
iy = 21
maxTr = 3.1
minTr = 0.1
maxTn = maxTr
minTn = minTr
Tr =np.linspace(minTr, maxTr, num=ix, endpoint=True)
tr_t = Tr*2  ## tr must be multiple of dt
Tn = np.linspace(minTn, maxTn, num=iy, endpoint=True)
tn = 3  # tn = 3 hours
K = tn/n_sub

#Q_cso_SI = 8.0 # mm/h 
Q_cso_US = 0.3 # in/h (US unit) it is normalized by area
cutoff = Q_cso_US*all_area/43560.0  #
cso_xy= np.zeros((ix,iy))  # Mgal
peak= np.zeros((ix,iy))    # cfs
F_xy= np.zeros((ix,iy))
#%%
# convert the reusults into SI unit
max_CSO_V = []
max_CSO_perc = []
max_Peak = []
max_Peak_perc = []
for location in range(0,(n_sub+1)):
    for i in range(0,ix):
        for m in range(0,iy):
            rain_i = Tn[m]*DR/(tn*bern_h)# i increases with tr    
            tr = Tr[i]*bern_h/(DR*rain_i)  # tr must be a multiple of dt; i.e. tr/dt must be integer
            storm_v = rain_i*tr*25.4 #mm
            rain= list(np.ones(int(round(tr/dt)))*rain_i)+list(np.zeros(int(20/dt)))                
            End_T= len(rain)*dt
            inflow=np.zeros(len(rain))
            Q_t,Q_all,cso_m, F = cso_Mgal(tr,K,x,cutoff,inflow,rain,all_area,tot_area,imperv_perc,gi_s_p,DR,Sf,L,dt,End_T,n_sub,location, manning_n,bern_h)
#            cso_m: Million gal
            cso_xy[i,m]=cso_m*1e6*0.0037854*0.1/storm_v   # mm
            Q_t = Q_t * 0.028316 # converting cfs to cms
#            Q_all = Q_all * 0.028316 # cms

            peak[i,m]=max(Q_t)*0.1*3600/(rain_i*25.4)  #mm/mm
            F_xy[i,m]=F*1e6*0.0037854*0.1/storm_v   # mm
#            print "K= %s, tr= %s Peak Flow = %s" %(k,tr,peak[i,m])
    ## 3D Plot (Tr - Tn)
    plt.rcParams.update({'font.size': 22})
    fig = plt.figure(figsize=(10, 6))
    
#    plt.figure(figsize=(8, 4))  # Tr : Tn
#    ax = plt.gca()    
#    ax.set_yscale('log')
#    ax.set_xscale('log')
    
    
    X, Y = np.meshgrid(Tr, Tn)
    zs=cso_xy
#    print "Max CSO  %s m3/ha" %max(max(cso_xy.tolist()))

    zs=np.array(map(list, zip(*zs)))
    Zc = zs.reshape(X.shape)
    norm = colors.Normalize(vmin = 0, 
                            vmax = 0.5)
    v = np.linspace(0,0.5, num=6, endpoint=True)    
#    im = plt.imshow((Zc),cmap='Blues',extent=(minTr, maxTr, minTn, maxTn), origin ='lower') # drawing the function

    im = plt.imshow((Zc),cmap='Blues',extent=(minTr, maxTr, minTn, maxTn), origin ='lower', norm = norm) # drawing the function
    plt.colorbar(im, ticks=v) # adding the colobar on the right
    plt.xlabel('Tr')
    plt.ylabel('Tn')
    ticks = [0.1]
    ticks.extend(list(np.arange(1,maxTr,step=1)))
    plt.xticks(ticks)
    plt.yticks(ticks)        
#        plt.title('CSO_Reduction (CSO(mm)')        
    plt.tight_layout()
    plt.grid(True)

    if location == 0:
        cso_no=Zc
        plt.savefig("3D_CSO_No_GI_tgi_2")
        plt.close(fig)  

    else:
        plt.savefig("3D_CSO_Location%s_tgi_2" %location)
        plt.close(fig)  
        
#### CSO Reduction(mm) plot  
        fig = plt.figure(figsize=(10, 6))
        norm = colors.Normalize(vmin = 0, 
                                vmax = 0.1)
        v = np.linspace(0,0.1, num=6, endpoint=True)    
        Zcp = (cso_no-Zc)
#        im = plt.imshow(Zcp,cmap='RdYlBu',extent=(minTr, maxTr, minTn, maxTn), origin ='lower') # drawing the function

        im = plt.imshow(Zcp,cmap='RdYlBu',extent=(minTr, maxTr, minTn, maxTn), origin ='lower', norm = norm) # drawing the function
        plt.grid(True)
        plt.colorbar(im, ticks=v) # adding the colobar on the right
        plt.xlabel('Tr')
        plt.ylabel('Tn')
        ticks = [0.1]
        ticks.extend(list(np.arange(1,maxTr,step=1)))
        plt.xticks(ticks)
        plt.yticks(ticks)        
#        plt.title('CSO_Reduction (CSO(mm)')        
        plt.tight_layout()
        plt.grid(True)
        plt.savefig("3D_CSO_Reduction_vol_Location%s_tgi_2" %location)
        plt.close(fig)  
#        print "Max CSO reduction %s m3/h" %max(max(Zcp.tolist()))

#        temp = (cso_no - Zc)
#        max(temp.tolist())
#### CSO Reduction(%) plot  
        fig = plt.figure(figsize=(10, 6))
        norm = colors.Normalize(vmin = 0, 
                                vmax = 100)
        Z1 = (cso_no-Zc)/cso_no*100
        Z2 = np.nan_to_num(Z1)
#        print "Max CSO reduction %s m3/h" %max(max(Z2.tolist()))

        v = np.linspace(0,100, num=5, endpoint=True)    
        im = plt.imshow((Z2),cmap='RdYlBu',extent=(minTr, maxTr, minTn, maxTn), origin ='lower', norm = norm) # drawing the function
        plt.grid(True)
        plt.colorbar(im, ticks=v) # adding the colobar on the right
        plt.xlabel('Tr')
        plt.ylabel('Tn')
        ticks = [0.1]
        ticks.extend(list(np.arange(1,maxTr,step=1)))
        plt.xticks(ticks)
        plt.yticks(ticks)        
#        plt.title('CSO_Reduction (%)')        
        plt.tight_layout()
        plt.savefig("3D_CSO_Reduction_perc_Location%s_tgi_2" %location)
        plt.close(fig)  
        
        max_CSO_V.append(max(max(Zcp.tolist())))        
        max_CSO_perc.append(max(max(Z2.tolist())))        
    
## 3D Plot (Peak): peak discharge should be close to rainfall intensity i:12.7 mm/h
    fig = plt.figure(figsize=(10, 6))
    zs=peak # mm/mm
#    print "Max Peak Flow %s mm/h" %max(max(zs.tolist()))
    zs=np.array(map(list, zip(*zs)))
    Zp = zs.reshape(X.shape)
    # MAIN IDEA: Added normalisation using nanmin and nanmax functions
    norm = colors.Normalize(vmin = 0, 
                            vmax = 1)
    v = np.linspace(0, 1, 6, endpoint=True)            
#    im = plt.imshow((Zp),cmap='Greens',extent=(minTr, maxTr, minTn, maxTn), origin ='lower') # drawing the function
    im = plt.imshow((Zp),cmap='Greens',extent=(minTr, maxTr, minTn, maxTn), origin ='lower', norm = norm) # drawing the function
    
    plt.grid(True)
    plt.xlabel('Tr')
    plt.ylabel('Tn')
    plt.colorbar(im, ticks=v) # adding the colobar on the right
    ticks = [0.1]
    ticks.extend(list(np.arange(1,maxTr,step=1)))
    plt.xticks(ticks)
    plt.yticks(ticks)        
#        plt.tight_layout()
#        plt.title('Peak Flow (mm/h)')        
    plt.tight_layout()

    if location == 0:
        peak_no=Zp  # cfs
        plt.savefig("3D_Peak_No_GI_tgi_2")
        plt.close(fig)  
    else:

        plt.savefig("3D_Peak_Location%s_tgi_2" %location)
        plt.close(fig)  
                             
#### Peak Flow reduction Plot
        fig = plt.figure(figsize=(10, 6))
        norm = colors.Normalize(vmin = 0, 
                                vmax = 0.2)
        v = np.linspace(0, 0.2, 6 , endpoint=True)    
        Zd = (peak_no-Zp)        
        im = plt.imshow((Zd),cmap='RdYlGn',extent=(minTr, maxTr, minTn, maxTn), origin ='lower', norm =norm) # drawing the function
        plt.grid(True)
        plt.xlabel('Tr')
        plt.ylabel('Tn')
        plt.colorbar(im, ticks=v) # adding the colobar on the right
        ticks = [0.1]
        ticks.extend(list(np.arange(1,maxTr,step=1)))
        plt.xticks(ticks)
        plt.yticks(ticks)        
        plt.tight_layout()
#        print "Max Peak reduction %s mm/h" %max(max(Zd.tolist()))
#        plt.title('Peak Flow Reduction(mm/h)')     
        plt.savefig("3D_Peak_Reduction_flow_Location%s_tgi_2" %location)   
        plt.close(fig)  

#### Peak Flow reduction (%) Plot
        fig = plt.figure(figsize=(10, 6))
        norm = colors.Normalize(vmin = 1, 
                                vmax = 30)
        v = np.linspace(0, 30, 7 , endpoint=True)    
        Z1 = (peak_no-Zp)/peak_no*100
        Z2 = np.nan_to_num(Z1)
        im = plt.imshow(Z2,cmap='RdYlGn',extent=(minTr, maxTr, minTn, maxTn), origin ='lower', norm = norm) # drawing the function
        plt.grid(True)
        plt.xlabel('Tr')
        plt.ylabel('Tn')
        plt.colorbar(im, ticks=v) # adding the colobar on the right
        ticks = [0.1]
        ticks.extend(list(np.arange(1,maxTr,step=1)))
#        plt.title('Peak Flow Reduction(%)')

        plt.xticks(ticks)
        plt.yticks(ticks)        
        plt.tight_layout()
        plt.savefig("3D_Peak_Reduction_perc_Location%s_tgi_2" %location)   
        plt.close(fig)  
        
        max_Peak.append(max(max(Zd.tolist())))        
        max_Peak_perc.append(max(max(Z2.tolist())))    
    GN = tn*Q_cso_US*rain_i*DR/bern_h
    print('GN= ', GN)
#%% Contour
#        levels = np.arange(0, 15, 0.5)

## Contour Map
#        CS = plt.contour(X, Y, Z2, levels, cmap='RdYlGn', norm=norm, linewidths=2, extent=(.1, 2, 2, 0.1))
#        manual_locations = [(-1, -1.4), (-0.62, -0.7), (-2, 0.5), (1.7, 1.2), (2.0, 1.4), (2.4, 1.7)]
#        ax.clabel(CS, inline=1, fontsize=10, manual=manual_locations)
# Thicken the zero contour.
#        zc = CS.collections[6]
#        plt.setp(zc, linewidth=4)
#        plt.clabel(CS, levels[2::4],  # label every second level
#                  inline=1, fmt='%1.1f', fontsize=14)
#        plt.xticks(np.arange(0.2, 2.1, step=0.2))
#        plt.yticks(np.arange(0.2, 2.1, step=0.2))        
#
#        # latex fashion title
#        plt.title('Peak Flow Reduction(%)')
##        fig.colorbar(CS, shrink=0.5, aspect=4, ticks=v)         
#
#        fig.tight_layout()
#        plt.show()
#


#%%  Color Map
        
        
        



