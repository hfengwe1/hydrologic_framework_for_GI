# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 17:45:42 2017

@author: fengwei
"""
from __future__ import division
import math as ma

def f(d,alpha):
    return alpha*d**(5/3)

#dy = RK4(lambda t, y:alpha*y**(5/3))
#t,y,dt = 0., 0.2, .01
#
#while t <= 10:
#    if abs(round(t) - t )  < 1e-5:
#        print("y(%2.1f)\t = %4.6f " % (t, y))
#    t, y = t+dt, y-dy(t,y,dt)  
    
def RK4(f):
    return lambda t,y,dt:(
            lambda dy1:(
            lambda dy2:(
            lambda dy3:(
            lambda dy4:( dy1 + 2*dy2 + 2*dy3 + dy4 )/6
            )( dt * f(t+dt,   y+dy3  ) )
            )( dt * f(t+dt/2, y+dy2/2) )
            )( dt * f(t+dt/2, y+dy1/2) )
            )( dt * f(t, y           ) )
            
#from math import sqrt
#def theory(t): return (t**2 + 4)**2 /16
#
#dy = RK4(lambda t, y:t*sqrt(y))
#
#t,y,dt = 0., 1., .1
#while t <= 10:
#    if abs(round(t) - t )  < 1e-5:
#        print("y(%2.1f)\t = %4.6f \t error: % 4.6g" % (t, y, abs(y- theory(t))))
#    t, y = t+dt, y+dy(t,y,dt)