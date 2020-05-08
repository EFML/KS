#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 13:59:39 2017

@author: justinpringle
"""

import KS
import numpy as np
import matplotlib.animation as anim
import matplotlib.pyplot as plt
import datetime as dat
import time

def init():
    line.set_data([],[])
    return line,

def animate(i):
    
#    ax.cla()
    U = uStore[i]
    V = vStore[i]
    fac = int(np.floor_divide(len(xp),nt))
    sp = np.sqrt(U**2+V**2)
    
#    q = ax.quiver(x,y,U,V,scale=400)
    q.set_UVC(U,V)
#    s = ax.contourf(X,Y,sp,alpha=0.5)
#    print(xp[i],yp[i])
    ax.set_title('num %d'%i)
    line.set_data(xp[:i*fac],yp[:i*fac])
    line2.set_data(xp2[:i*fac],yp2[:i*fac])
    line3.set_data(xp3[:i*fac],yp3[:i*fac])
    line4.set_data(xp4[:i*fac],yp4[:i*fac])
    
    return line, q, line2,line3,line4,

def getVel(xp,yp,U,V,x,y):
    
    dx = x[0,1]-x[0,0]
    dy = y[1,0]-y[0,0]
    
    indX = int(np.floor_divide(xp,dx))
    ddx = np.remainder(xp,dx)/dx
    indY = int(np.floor_divide(yp,dy))
    ddy = np.remainder(yp,dy)/dy
#    print(ddy)
    u11 = U[indY,indX]
    u12 = U[indY+1,indX]
    u22 = U[indY+1,indX+1]
    u21 = U[indY,indX+1]
    
    u = ddy*(1-ddx)*u12 + ddx*ddy*u22 + ddx*(1-ddy)*u21 + (1-ddx)*(1-ddy)*u11
    
    v11 = V[indX,indY]
    v12 = V[indX,indY+1]
    v22 = V[indX+1,indY+1]
    v21 = V[indX+1,indY]
    
    v = ddy*(1-ddx)*v12 + ddx*ddy*v22 + ddx*(1-ddy)*v21 + (1-ddx)*(1-ddy)*v11
    
    return u,v

def updateParticle(t,dt,x,y,U,V,xp,yp):
    
    dtInt = 0.05 #integration time step
    tr = dtInt
    
    tempX = [xp[-1]]
    tempY = [yp[-1]]
    while tr < dt:
        tprime = t+0.5*dtInt
        u,v = getVel(tempX[-1],tempY[-1],U,V,x,y)
        xd = tempX[-1] + 0.5*dtInt*u
        yd = tempY[-1] + 0.5*dtInt*v
        if xd>x[0,-1] or yd >y[-1,0]:
            break
        elif xd<x[0,0] or yd <y[0,0]:
            break
        uu,vv = getVel(xd,yd,U,V,x,y)
    
        dx = uu*dtInt
        dy = vv*dtInt
        if tempX[-1]+dx>x[0,-1] or tempY[-1]+dy >y[-1,0]:
            break
        elif tempX[-1]+dx<x[0,0] or tempY[-1]+dy <y[0,0]:
            break
        tempX.append(tempX[-1]+dx)
        tempY.append(tempY[-1]+dy)
        tr+=dtInt
    
    xp.extend(tempX)
    yp.extend(tempY)
#    xp.append(xp[-1]+dx)
#    yp.append(yp[-1]+dy)
    
    return xp,yp
if __name__ == '__main__':
    
    #time in seconds
    dt = 55
    t=0
    nt = int(1.15e4/dt)
    
    #init
    dx = 780
    dy = 780
    nx = 100
    ny = 100
    Ms = 50
    Ml = 50
    L1 = dx*100
    Lg = dx
    Ln = 80
    #grid
    x,y = KS.genGrid(dx,dy,nx,ny)
    
    vel = KS.velocity(x,y,Ms,Ml,L1,Lg,Ln,t)
    vel.energySpec(0.2,0,250)
    
    uStore = []
    vStore = []
    
    #test particle tracking
    col,row = np.random.randint(0,high=x.shape[1]),np.random.randint(0,high=y.shape[0])
    xp = [x[0,col-1]]
    yp = [y[row,0]]
    
    xp2 = [x[0,col]]
    yp2 = [y[row,0]]
    xp3 = [x[0,col-2]]
    yp3 = [y[row,0]]
    xp4 = [x[0,col-3]]
    yp4 = [y[row,0]]
    
    now = time.mktime(dat.datetime.now().timetuple())
    
#    u,v = vel.calcLargeScale()
#    us,vs = vel.calcSmallScale()
    for i in range(nt):
        u,v = vel.calcLargeScale()
        us,vs = 0,0#vel.calcSmallScale()
        U = u+us
        V = v+vs
        uStore.append(U)
        vStore.append(V)
        
#        vel.updateTime(dt)
        
        xp,yp = updateParticle(i,dt,x,y,U,V,xp,yp)
        xp2,yp2 = updateParticle(i,dt,x,y,U,V,xp2,yp2)
        xp3,yp3 = updateParticle(i,dt,x,y,U,V,xp3,yp3)
        xp4,yp4 = updateParticle(i,dt,x,y,U,V,xp4,yp4)
        
    end = time.mktime(dat.datetime.now().timetuple())
    print('time: %.2f'%(end-now))
    X,Y = np.meshgrid(x,y)
    
    fig,ax = plt.subplots()
    ax.set_xlim(np.min(xp)-10*dx,np.max(xp)+10*dx)
    ax.set_ylim(np.min(yp)-10*dy,np.max(yp)+10*dy)
#    line, = ax.plot(xp,yp,linewidth=5)
    q = ax.quiver(x,y,uStore[0],vStore[0],scale=400,alpha=0.5)
#    s = ax.contourf(X,Y,np.zeros(X.shape),alpha=0.0)
    line, = ax.plot(xp,yp)
    line2, = ax.plot(xp2,yp2)
    line3, = ax.plot(xp3,yp3)
    line4, = ax.plot(xp4,yp4)
    
    aniM = anim.FuncAnimation(fig,animate,frames=nt,blit=False,init_func=init)
    plt.show()
        
        
        
        
        