#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 04:43:09 2017

@author: justinpringle

cloud dispersion model
"""

import numpy as np

class particle():
    '''
    particle class.
    each particle has an x and y coordinate this will then be interpolated 
    onto computational grid
    '''
    
    def __init__(self,x,y):
        self.x=x
        self.y=y
        
    def getLoc(self,X,Y):
        '''
        takes grid and interpolates position of particle.
        also used in biliner interpolation scheme
        '''
        dx = x[0,1]-x[0,0]
        dy = y[1,0]-y[0,0]
    
        indX = int(np.floor_divide(xp,dx))
        ddx = np.remainder(xp,dx)
        indY = int(np.floor_divide(yp,dy))
        ddy = np.remainder(yp,dy)
        
class cdm():
    '''
    cloud velocity has three components:
        a) U (mean velocity).
        b) u(x,t) - turbulent velocity field
        c) ur - random component modelling affects above wave numbers > kG
    Trajectory is calculated using Lagrangian integration.
    KS1 -> small scale turbulence, determineshow the cloud grows.
    Once cloud has grown to Lg it is split into several smaller clouds.
    This can cause memory issues so whole distribution of clouds is re-discretized
    involving splitting some clouds but merging others.
    '''
    
    def __init__(self,x,y,Lr,numPart):
        '''
        keep track of three things for cloud
        
        '''
        
        self.x = x
        self.y = y
        self.amount = amount
        self.centroid = centroid
        self.radius = Lr
        self.numPart = numPart
        
#    def arrangeParticles(self,x,y):
