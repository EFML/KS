#/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 05:30:42 2017

@author: justinpringle

Kinematic Simulations... The BEGINNING...

This script is used to plot the velocity function:

    u(x,t) = sum_(n=1)^Nk[An*cos(kn*x +wn*t)+Bn*sin(kn*x+wn*t)]
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import scipy.interpolate as interp

def genGrid(dx,dy,nx,ny):
    '''
    generate some grid with spacing dx,dy
    '''
    x = np.linspace(0,nx*dx,nx)
    y = np.linspace(0,ny*dy,ny)
    
    xg,yg = np.meshgrid(x,y,sparse=True)
    Y,X = np.mgrid[yg[0,0]:xg[0,-1]:100j,yg[0,0]:yg[-1,0]:100j]
    
    return X,Y

def genWaveNumber(Nk,L1,L2,first=False):
    '''
    generates array of wave numbers = to the number of modes.
    for large scale:
        2pi/Lg<k<2pi/Ls (length of grid and length scale)
    for small scale:
        2pi/Lr<k<2pi/Lg (length od cloud radius)
    k is geometric:
        kn = k1*alpha^(n-1)
        alpha = (L/neta)^(1/(nk-1))
    '''
    #initialize
    eta = L2
    k_0 = 2*np.pi/L1

    k = np.zeros(Nk)
    k[0] = k_0
    for i in range(1,Nk):
        alp = (L1/eta)**(1.0/(Nk))
        k[i] = k_0*(alp)**(i)
##    print(k)
    return k

def energySpec(Ls,U1,V1,k):
    '''
    uses argument in Fung 08
    L is the integral length scale, U1, V1 are vel components to calc energy.
    k are the wave numbers
    returns interpolated function
    '''
    
    u = np.sqrt(U1**2+V1**2)
    
    g1 = np.pi**2*(u**2)/(4*Ls**2)
    g2 = np.pi/(4*Ls)
    
#    k1 = 2*np.pi/L1
#    kn = 2*np.pi/Ln
#    k = np.zeros(100)
#    k[0] = k1
#    
#    for i in range(1,Nk):
#        alp = (L1/L2)**(1.0/Nk)
#        k[i] = k1*alp**i
#    #k = np.linspace(k1,kn,100)
#    E = np.zeros(k.shape[0])
    
    A1 = g1*k**3
    B1 = (g2**2+k**2)**3
    E = A1/B1
#    for i,j in enumerate(k):
#        print j,E[i]
#    for i,e in enumerate(E):
#        print(e,k[i])
#    eFunc = interp.interp1d(k,E)    
#    kFunc = interp.interp1d(E,k)
    
    return E

def getEnergy(k,eFunc,kFunc,test=True):
    '''
    test assumes energy is gausian distributed
    '''
    
    if test:
        E = 10*st.norm.pdf(k,loc=1,scale=0.5)
    else:
        E = eFunc(k)/100
        
    return E

def sumVel(Nk,x,y,t,L1,L2,U1,V1):
    '''
    sum velocity components
    Nk is the number of modes
    x,y are grid points
    kn,An,Bn are vectors
    '''
    #select k values to calc deltaK etc.
#    L1 = 10**5#
#    L2=4*np.max([x[0,1]-x[0,0],y[1,0]-y[0,0]])
#    Lr=80
#    print(L2)
    k = genWaveNumber(Nk,L1,L2)
    en = energySpec(L,U1,V1,k)
    
    lamm = 0.3
    #loop over modes and sum values
    u = np.zeros((y.shape[0],x.shape[1]))
    v = np.zeros((y.shape[0],x.shape[1]))
    
    for i in range(Nk):
        #randomly select phi
        phiN = np.deg2rad(np.random.uniform(0,360))
        
        #calc delta k
        if i == 0:           
            deltK = (k[1]-k[0])/2
        elif i == Nk-1:            
            deltK = (k[Nk-1]-k[Nk-2])/2
        else:
            deltK = (k[i+1]-k[i-1])/2
            
        #get energy given wave number
        kn = k[i]
        E = en[i]#getEnergy(kn,eFunc,kFunc,test=False)
#        print(E)
        if E == 0:
            continue
               
        stdAB = np.sqrt(E*deltK)
        
        An = st.norm.cdf(np.random.uniform(),scale=stdAB)
        Bn = st.norm.cdf(np.random.uniform(),scale=stdAB)
        
        
        wn = lamm*np.sqrt(kn**3*E)
        
        #dot product between space vector and wave number vector
        knDx = kn*np.sin(phiN)*x+kn*np.cos(phiN)*y
        u+=(An*np.cos(phiN)*np.cos(knDx+wn*t) + Bn*(-np.cos(phiN))*np.sin(knDx+wn*t))
        v+=(An*(-np.sin(phiN))*np.cos(knDx+wn*t) + Bn*(np.sin(phiN))*np.sin(knDx+wn*t))
        
    return u,v

if __name__ == '__main__':
    
    #first generate the grid
    dx = 780    
    dy = 780
    nx = 100
    ny = 100
        
    grid_x,grid_y = genGrid(dx,dy,nx,ny)
    
    #energy spectrum
    #   assume integral length scale of 250m as Fung 08, U=1m/s
    U1 = 1
    V1=0
    L = 250
    
    #now generate some velocity field
    Nk = 50
    nLS = 10
    L1 = 10**5#
    L2=np.max([grid_x[0,1]-grid_x[0,0],grid_y[1,0]-grid_y[0,0]])
    u1,v1 = sumVel(Nk,grid_x,grid_y,0,L1,L2,U1,V1)
    
    #small scale
    Lr=80
    u2,v2 = sumVel(Nk,grid_x,grid_y,0,L2,Lr,U1,V1)
    
    u = u1+u2
    v = v1+v2
    
    Y,X = np.mgrid[grid_x[0,0]:grid_x[0,-1]:100j,grid_y[0,0]:grid_y[-1,0]:100j]
    sp = np.sqrt(u*u+v*v)
#    plt.quiver(grid_x[:,0:-1:3],grid_y[0:-1:3,:],u[0:-1:3,0:-1:3],v[0:-1:3,0:-1:3],scale=100)
#    plt.quiver(grid_x,grid_y,u,v,scale=400)
#    plt.contourf(X,Y,sp,alpha=0.5)
#    plt.streamplot(grid_x,grid_y,u,v,linewidth=sp/10,density=2)
#    plt.streamplot(grid_x,grid_y,u,v,linewidth=sp)
#    plt.show()