# -*- coding: utf-8 -*-
"""
Created on Fri Apr 07 09:43:26 2017
Gravity Sewer - no pumping 
@author: fengwei
"""
from __future__ import division
def Manning(n,Sf,Rx,Ax):  
    Q=1.49/n*Sf**0.5*Rx**(2/3)*Ax
    return Q
