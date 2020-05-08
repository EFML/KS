#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 08:57:52 2017

@author: justinpringle

KS class scripts

This may be an easier way of keeping tracks of specific sims i.t.o temporal 
structure etc.
"""

import numpy as np
import scipy.stats as st

def genGrid(dx,dy,nx,ny):
    '''
    generate some grid with spacing dx,dy
    '''
    x = np.linspace(0,nx*dx,nx)
    y = np.linspace(0,ny*dy,ny)
    
    xg,yg = np.meshgrid(x,y,sparse=True)
    #Y,X = np.mgrid[yg[0,0]:xg[0,-1]:100j,yg[0,0]:yg[-1,0]:100j]
    
    return xg,yg

class velocity():
    '''
    constructs spatial velocity structures from a number of modes
    this class calculates large scale and small scale vels
    '''
    def __init__(self,x,y,Ms,Ml,L1,Lg,Ln,t):
        '''
        takes in grid (x,y)
        Number of modes -> Ms = small scale, Ml = large scale
        L1 -> starting wave length (large scale)
        Lg -> grid scale
        Ln -> end wave length (small scale)        
        '''
        #store some stats as lists in tuple smaller scale first
        self.An = ([],[])
        self.Bn = ([],[])
        self.phiN =([],[])
        self.kn = ([],[])
        self.wN = ([],[])
        self.Nk = (Ms,Ml)
        self.L1 = L1
        self.Lg = Lg
        self.Ln = Ln
        self.t=t
        self.x=x
        self.y=y
        
    def updateTime(self,dt):
        '''
        
        '''
        self.t+=dt
    
    def energySpec(self,u1,u2,Ls):
        '''
        calculates the default homogenous isotropic energy spectrum
        '''
        u = np.sqrt(u1**2+u2**2)
        self.g1 = np.pi**2*(u**2)/(4*Ls**2)
        self.g2 = np.pi/(4*Ls)
    
    def energyFunc(self,k):
        
        A1 = self.g1*k**3
        B1 = (self.g2**2+k**2)**3
        E = A1/B1
        
        return E
    
    def calcVelocity(self,N,alp,k0,large=True):
        '''
        grunt work
        '''
        lamm = 30
        
        u = np.zeros((self.y.shape[0],self.x.shape[1]))
        v = np.zeros((self.y.shape[0],self.x.shape[1]))
        
        if self.t==0:
            #then the beginning
            for i in range(N):
                phiN = np.deg2rad(np.random.uniform(0,360))
                k = k0*alp**i
                #calc delta k
                if i == 0:
                    k1 = k0*alp**(i+1)
                    deltK = (k1-k0)/2
                elif i == N-1:
                    kN = k0*alp**i
                    kN2 = k0*alp**(i-1)
                    deltK = (kN-kN2)/2
                else:
                    k1 = k0*alp**(i+1)
                    k2 = k0*alp**(i-1)
                    deltK = (k1-k2)/2
                
                E = self.energyFunc(k)
                
                stdAB = np.sqrt(E*deltK)
                An = st.norm.cdf(np.random.uniform(),scale=stdAB)
                Bn = st.norm.cdf(np.random.uniform(),scale=stdAB)
                
                wn = lamm*np.sqrt(k**3*E)
                
                #dot product between space vector and wave number vector
                knDx = k*np.sin(phiN)*self.x+k*np.cos(phiN)*self.y
                
                u+=(An*np.cos(phiN)*np.cos(knDx+wn*self.t) + Bn*(-np.cos(phiN))*np.sin(knDx+wn*self.t))
                v+=(An*(-np.sin(phiN))*np.cos(knDx+wn*self.t) + Bn*(np.sin(phiN))*np.sin(knDx+wn*self.t))
                
                if large:
                    self.An[1].append(An)
                    self.Bn[1].append(Bn)
                    self.phiN[1].append(phiN)
                else:
                    self.An[0].append(An)
                    self.Bn[0].append(Bn)
                    self.phiN[0].append(phiN)
#                self.An[1].append(An)
#                self.Bn[1].append(Bn)
#                self.phiN[1].append(phiN)
        else:
            for i in range(N):
                phiN = self.phiN[1][i]
                k = k0*alp**i
                #calc delta k
                if i == 0:
                    k1 = k0*alp**(i+1)
                    deltK = (k1-k0)/2
                elif i == N-1:
                    kN = k0*alp**i
                    kN2 = k0*alp**(i-1)
                    deltK = (kN-kN2)/2
                else:
                    k1 = k0*alp**(i+1)
                    k2 = k0*alp**(i-1)
                    deltK = (k1-k2)/2
                
                E = self.energyFunc(k)
                stdAB = np.sqrt(E*deltK)
                if large:
                    An = self.An[1][i]
                    Bn = self.Bn[1][i]
                else:
                    An = self.An[0][i]
                    Bn = self.Bn[0][i]
                
                wn = lamm*np.sqrt(k**3*E)
                
                #dot product between space vector and wave number vector
                knDx = k*np.sin(phiN)*self.x+k*np.cos(phiN)*self.y
                u+=(An*np.cos(phiN)*np.cos(knDx+wn*self.t) + Bn*(-np.cos(phiN))*np.sin(knDx+wn*self.t))
                v+=(An*(-np.sin(phiN))*np.cos(knDx+wn*self.t) + Bn*(np.sin(phiN))*np.sin(knDx+wn*self.t))
            
                
        return u,v            
                
    def calcLargeScale(self):
        '''
        calcs the large scale velocity field
        '''
        #wave numbers and energy
        N = self.Nk[1]
        alp = (float(self.L1)/self.Lg)**(1.0/N)
        k0 = 2*np.pi/self.L1
        
        u,v = self.calcVelocity(N,alp,k0)
        
        return u,v
        
                
        
    def calcSmallScale(self):
        '''
        calcs the small scale velocity field
        '''
        #wave numbers and energy
        N = self.Nk[0]
        
        alp = (float(self.Lg)/self.Ln)**(1.0/N)
        
        k0 = 2*np.pi/self.Lg
        
        u,v = self.calcVelocity(N,alp,k0,large=False)
        
        return u,v
        
        
        
        
        
        
        
        
        