#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:36:05 2023

@author: charlesarnold
"""

print('Homework #1, due Sunday Jan 20')
import numpy as np
from numpy import pi,sin,exp,cos
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy.signal import find_peaks
from scipy import signal
from PIL import Image

# In[9]:



def integral(P,a):
    t0 = 0
    t1 = 10
    dt = 0.01
    d = 1
    B = 2
    K=2
    m = int((t1-t0)/dt)
    time = np.zeros(m+1)
    P_0 = np.zeros(m+1)
    time[0]= t1
    P_0[0] = P
    
    for i in range(1, m+1):
       t = time[i-1]
       P_p = P_0[i-1]
       P_dot = a+(B*P_p**(100))/(K**(100)+P_p**(100))-d*P_p
       time[i] = t + dt
       P_0[i] = P_p+P_dot*dt
      
    return time, P_0

a = np.linspace(0,4,6)

c = np.linspace(1,5,6)
for j in range(len(c)):
    plt.figure()
    for i in range(len(a)):
        plt.plot(*integral(round(c[j]),a[i]), linestyle = "-", label = f'P = {round(c[j],3)}; a = {round(a[i])}') 
        plt.xlabel("Time") 
        plt.ylabel("P") 
        plt.legend()
        plt.title("chaotic coding")  

# In[9]:
# In[9]:
r1 = np.linspace(-1,2,1000)
y = []
for i in range(len(r1)):
    y.append(integral(4,r1[i])[1][-1])
plt.plot(r1, y)
plt.title('a vs xf')
plt.xlabel('a')
plt.ylabel('xf')
# In[9]:

# In[9]:

# In[9]:
Ex = []; Et = []; Edt = 0.009;
Rkx = []
Rkt = []
Rky = []
Rkdt = 0.009
tf = 10
rt = 5

tt=0;xt=0.1;
while(tt<tf):
    vt = rt*xt-xt/(1+xt**2) 
    xt = xt+Edt*vt
    tt = tt+Edt
    Ex.append(xt)
    Et.append(tt)


    

    


# In[9]:
tt=0;xt=0.1;
while(tt<tf):
    x1 = xt
    k1 = (rt*x1)
    x2 = xt+Rkdt*k1/2
    k2 = (rt*x2)
    x3 = xt+Rkdt*k2/2
    k3 = (rt*x3)
    x4 = (xt+Rkdt*k3/2)
    k4 = (rt*x4)
    
    vt = (k1+2*k2+2*k3+k4)/6   
    xt = xt+Rkdt*vt
    tt = tt+Rkdt
    Rkx.append(xt)
    Rkt.append(tt)

# In[9]:
tt=0;xt=0.1;
while(tt<tf):
    x1 = xt
    k1 = (rt*x1-x1/(1+x1**2))
    x2 = xt+Rkdt*k1/2
    k2 = (rt*x2-x2/(1+x2**2))
    x3 = xt+Rkdt*k2/2
    k3 = (rt*x3-x3/(1+x3**2))
    x4 = xt+Rkdt*k3/2
    k4 = (rt*x4-x4/(1+x4**2))
    
    vt = (k1+2*k2+2*k3+k4)/6
    xt = xt+Rkdt*vt
    tt = tt+Rkdt
    Rkx.append(xt)
    Rkt.append(tt)

# In[9]:
plt.plot(Et,Ex, color = 'red', label = f'Edt = {Edt}')
plt.xlabel("time (s)")
plt.ylabel("X(t) (m)")
plt.legend()
plt.figure()
plt.plot(Rkt, Rkx, color = 'blue', label = f'Rkdt = {Rkdt}')
plt.xlabel("time (s)")
plt.ylabel("X(t) (m)")
plt.legend()

# In[9]:

# In[9]:

# In[9]:

# In[9]:# In[9]:
def integralE(N,r,k):
    t0 = 0
    t1 = 10
    dt = 0.01
    n = int((t1-t0)/dt)
    time = np.zeros(n+1)
    N_0 = np.zeros(n+1)
    time[0]= t1
    N_0[0] = N
    for i in range(1, n+1):
       t = time[i-1]
       N_n = N_0[i-1]
       N_dot = r*N_n*(1-N_n/k) 
       time[i] = t + dt
       N_0[i] = N_n+N_dot*dt
      
    return time, N_0

# In[9]:
tt=0;xt=0.1;
while(tt<tf):
    vt = rt*xt-xt/(1+xt**2) 
    xt = xt+Edt*vt
    tt = tt+Edt
    Ex.append(xt)
    Et.append(tt)

# In[9]:
    
Ex = []; Et = []; Edt = 0.009;
tt=0;xt=0.1;
while(tt<tf):
    vt = rt*xt-xt/(1+xt**2) 
    xt = xt+Edt*vt
    tt = tt+Edt
    Ex.append(xt)
    Et.append(tt)


# In[9]:





# In[9]:
def coupled_rk4(dxdt, dydt, x_0, y_0, dt, tf):
    
    def f(xy,t):
        x = xy[0]
        y = xy[1]
        return np.array([dxdt(x,y,t), dydt(x,y,t)], float)
    
    t_arr = np.arange(0, tf, dt)
    x_arr = []
    y_arr = []
    
    x_arr.append(x_0)
    y_arr.append(y_0)
    
    xy = np.array([x_0, y_0], float)
    
    for t in t_arr:
        k1 = f(xy, t) * dt
        k2 = f(xy + 0.5*k1, t + 0.5*dt) * dt
        k3 = f(xy + 0.5*k2, t + 0.5*dt) * dt
        k4 = f(xy + k3, t + dt) * dt
        xy += (k1 + 2*k2 + 2*k3 + k4)/6
        x_arr.append(xy[0])
        y_arr.append(xy[1])
    
    return np.array(x_arr), np.array(y_arr)
# In[9]:

# In[9]:
x=[]; y=[]; vx=[]; vy=[];
x0=-10;xf=10;dx=2;
y0=-10;yf=10;dy=2;



for xt in range(x0,xf,dx):
    for yt in range(y0,yf,dy):
        x.append(xt)
        y.append(yt)
        vx.append(yt)
        vy.append(-2*xt-3*yt)
        
plt.quiver(x,y,vx,vy, color='grey', width=0.003)

points = [(1,0), (0,1), (2,1), (-2,-2)]
x = np.arange(-10,10,0.2)
y = x
plt.plot(-x/2,y,label=r'$\vec{v}_1$',linestyle='-.')
plt.plot(-x,y,label=r'$\vec{v}_2$',linestyle='-.')
     
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')  
for point in points:
    line = plt.plot(*coupled_rk4(dxdt, dydt, point[0], point[1], dt, tf), label=f'$x_0 = {{{point[0]}}}$, $y_0 = {{{point[1]}}}$')
    plt.scatter(point[0], point[1], color=line[0].get_color())
plt.title(r'$\dot{x} = y$, $\dot{y} = -2x-3y$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.legend()
# In[9]:



#%% Problem 1

dt = 0.01
tf = 5

x, y = np.meshgrid(np.linspace(-3, 3, 15), np.linspace(-3, 3, 15))
t = np.arange(0, tf, dt)

def dxdt(x,y,t): return y

def dydt(x,y,t): return -2*x-3*y

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(1,0), (0,1), (2,1), (-2,-2)]

plt.figure(dpi=300)
#plt.style.use('default')
#plt.style.use('fivethirtyeight')
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
plt.quiver(x, y, X, Y, color='grey', width=0.003)
for point in points:
    line = plt.plot(*coupled_rk4(dxdt, dydt, point[0], point[1], dt, tf), label=f'$x_0 = {{{point[0]}}}$, $y_0 = {{{point[1]}}}$')
    plt.scatter(point[0], point[1], color=line[0].get_color())
x = np.arange(-10,10,0.1)
y = x
#plt.plot(-x/2,y,label=r'$\vec{v}_1$',linestyle='-.')
#plt.plot(-x,y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = y$, $\dot{y} = -2x-3y$', fontname='Batang')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.legend(loc='upper right')
#leg = plt.legend(frameon=True)
#leg.get_frame().set_edgecolor('black')
plt.show()
# In[9]:
    
# In[9]:
dt = 0.01
tf = 5

x, y = np.meshgrid(np.linspace(-3, 3, 15), np.linspace(-3, 3, 15))
t = np.arange(0, tf, dt)

def dxdt(x,y,t): return y

def dydt(x,y,t): return -2*x-3*y

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(1,0), (0,1), (2,1), (-2,-2)]

plt.figure(dpi=300)
#plt.style.use('default')
#plt.style.use('fivethirtyeight')
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
plt.quiver(x, y, X, Y, color='grey', width=0.003)
for point in points:
    line = plt.plot(*coupled_rk4(dxdt, dydt, point[0], point[1], dt, tf), label=f'$x_0 = {{{point[0]}}}$, $y_0 = {{{point[1]}}}$')
    plt.scatter(point[0], point[1], color=line[0].get_color())
x = np.arange(-10,10,0.1)
y = x
plt.plot(-x/2,y,label=r'$\vec{v}_1$',linestyle='-.')
plt.plot(-x,y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = y$, $\dot{y} = -2x-3y$', fontname='Batang')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.legend(loc='upper right')
#leg = plt.legend(frameon=True)
#leg.get_frame().set_edgecolor('black')
plt.show()
# In[9]:

# In[9]:

# In[9]: oscilator from -3,3
    
dt = 0.01
tf = 5

x, y = np.meshgrid(np.linspace(-3, 3, 15), np.linspace(-3, 3, 15))
t = np.arange(0, tf, dt)

r = 1

def dxdt(x,y,t): return y

def dydt(x,y,t): return r*y-x

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(1,0), (0,1), (2,1), (-2,-2)]

plt.figure(dpi=300)
#plt.style.use('default')
#plt.style.use('fivethirtyeight')
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
plt.quiver(x, y, X, Y, color='grey', width=0.003)
for point in points:
    line = plt.plot(*coupled_rk4(dxdt, dydt, point[0], point[1], dt, tf), label=f'$x_0 = {{{point[0]}}}$, $y_0 = {{{point[1]}}}$')
    plt.scatter(point[0], point[1], color=line[0].get_color())
x = np.arange(-10,10,0.1)
y = x
plt.plot(r*y,y,label=r'$\vec{v}_1$',linestyle='-.')
plt.plot(np.linspace(0,0,200),y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = y$, $\dot{y} = {r}*y-x$' f', r = {r}')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.legend(loc='upper right')
#leg = plt.legend(frameon=True)
#leg.get_frame().set_edgecolor('black')
plt.show()


# In[9]:
    
dt = 0.01
tf = 5

x, y = np.meshgrid(np.linspace(-5, 25, 30), np.linspace(-5, 25, 30))
t = np.arange(0, tf, dt)

r = -2

def dxdt(x,y,t): return y

def dydt(x,y,t): return r*y-x

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(-5,25), (-5,15), (-5,20), (-5,10)]

plt.figure(dpi=300)
#plt.style.use('default')
#plt.style.use('fivethirtyeight')
plt.plot([-5,25], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-5,25], color='black', linewidth=0.75, linestyle='--')
plt.quiver(x, y, X, Y, color='grey', width=0.003)
for point in points:
    line = plt.plot(*coupled_rk4(dxdt, dydt, point[0], point[1], dt, tf), label=f'$x_0 = {{{point[0]}}}$, $y_0 = {{{point[1]}}}$')
    plt.scatter(point[0], point[1], color=line[0].get_color())
x = np.arange(-10,10,0.1)
y = x
plt.plot(r*y,y,label=r'$\vec{v}_1$',linestyle='-.')
plt.plot(np.linspace(0,0,200),y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = y$, $\dot{y} = {r}*y-x$' f', r = {r}')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-6, 25)
plt.ylim(-1, 15)
plt.legend(loc='upper right')
#leg = plt.legend(frameon=True)
#leg.get_frame().set_edgecolor('black')
plt.show()

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]: