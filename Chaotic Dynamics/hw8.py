#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 00:15:34 2023

@author: charlesarnold
"""


# In[9]:
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

def n_coupled_rk4(diffeqs, initials, dt, tf):
    
    def f(r,t):
        return np.array([diffeq(*r,t) for diffeq in diffeqs], float)
    
    t_arr = np.arange(0,tf,dt)
    
    out_arr = []
    for initial in initials:
        out_arr.append([initial])
    
    r = np.array(initials, float)
    
    for t in t_arr:
        k1 = f(r, t)*dt
        k2 = f(r + 0.5*k1, t + 0.5*dt)*dt
        k3 = f(r + 0.5*k2, t + 0.5*dt)*dt
        k4 = f(r + k3, t + dt)*dt
        r += (k1 + 2*k2 + 2*k3 + k4)/6
        for i in range(len(r)):
            out_arr[i].append(r[i])
    
    return out_arr


def gen_points(start, end, n_points):
    points = []
    for i in np.linspace(start,end,n_points):
        points.append((i,0))
        points.append((0,i))
        points.append((i,i))
        points.append((i,-i))
    return points

def polar_quiver(r, theta, drdt, dthetadt):
    R = drdt(r,theta,t)
    Theta = dthetadt(r,theta,t)
    x, y = r*np.cos(theta), r*np.sin(theta)
    X, Y = np.cos(theta)*R-r*np.sin(theta)*Theta,np.sin(theta)*R+r*np.cos(theta)*Theta
    plt.quiver(x, y, X, Y, color='grey', width=0.003)

# In[9]:
# In[9]:number 1
# In[9]:
dt = 0.01
tf = 50

#question4

x, y = np.meshgrid(np.linspace(-5, 5, 13), np.linspace(-5, 5, 13))
t = np.arange(0, tf, dt)

r = 1 #1, 0,-1 show how the phase diagram changes 

def dxdt(x,y,t): return r*x+y-x**2

def dydt(x,y,t): return -x+r*y+2*x**2

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = gen_points(-3,3,4)

plt.figure(dpi=300)
#plt.style.use('default')
#plt.style.use('fivethirtyeight')
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
plt.quiver(x, y, X, Y, color='grey', width=0.003)
for point in points:
    line = plt.plot(*n_coupled_rk4([dxdt, dydt], [point[0], point[1]], dt, tf), color='tab:blue')
    plt.plot(*n_coupled_rk4([dxdt, dydt], [point[0], point[1]], -dt, -tf), color=line[0].get_color())
    
x = np.arange(-10,10,0.1)
y = x

plt.plot(x,-r*x+x**2,label=r'$\dot{x} = 0$', color = 'purple', linestyle='-.')
plt.plot(x,(x/r-2/r*x**2),label=r'$\dot{y} = 0$',linestyle='-.')
plt.title(r'$\dot{x} = r*x+y+x^2$, $\dot{y} = -x+r*y+2*x^2$ $, ' f' r = {r}$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-5, 5)
plt.ylim(-5, 5)


plt.show()
# In[9]:number 2 

dt = 0.01
tf = 25

#question4


r, theta = np.meshgrid(np.linspace(0, 4, 10), np.linspace(0, 2*np.pi, 50))
rr = 1
t = np.arange(0, tf, dt)



def drrdt(r,rr,theta,t): return (h**2 / (m**2 * r**3) - k) / m

def drdt(r,rr,theta,t): return rr

def dthetadt(r,rr,theta,t): return h / (m*r**2)



m = 1
h = 1
k = 1

points = [(0,0,0), (1,1,1), (2,2,2)]

R = drdt(r,rr,theta,t)
Theta = dthetadt(r,rr,theta,t)
x, y = r*np.cos(theta), r*np.sin(theta)
X, Y = np.cos(theta)*R-r*np.sin(theta)*Theta,np.sin(theta)*R+r*np.cos(theta)*Theta


plt.figure(dpi=300)
#plt.style.use('default')
#plt.style.use('fivethirtyeight')
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
plt.quiver(x, y, X, Y, color='grey', width=0.003)

for point in points:
    plt.scatter(point[0]*np.cos(point[2]), point[0]*np.sin(point[2]))
    r, rr, theta = n_coupled_rk4([drdt, drrdt, dthetadt], [point[0], point[1], point[2]], dt, tf)
    plt.plot(r*np.cos(theta), r*np.sin(theta), label=f'$r_0 = {{{point[0]}}}$')



#plt.plot(x,-1/(a*x**2)*(1-(b+1)*x),label=r'$\dot{x} = 0$',linestyle='-.')
#plt.plot(x,b/(a*x),label=r'$\dot{y} = 0$',linestyle='-.')
plt.title(r'$\ddot{r}=\frac{h^2}{m^2r^3}-k$, $\dot{\theta}=\frac{h}{mr^2}$'+' \n'+f'$m={{{m}}}$, $h={{{h}}}$, $k={{{k}}}$')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.legend()

plt.show()
# In[9]: number 3 

dt = 0.01
tf = 50

#question4

x, y = np.meshgrid(np.linspace(-5, 5, 15), np.linspace(-5, 5, 15))
t = np.arange(0, tf, dt)

a = 1
b = 1

def dxdt(x,y,t): return 1-(b+1)*x+a*x**2*y

def dydt(x,y,t): return b*x-a*x**2*y

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(1,1), (-3,-1), (-2.5,1), (1,-1), (0.1,0)]

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

plt.plot(x,1/(a*x**2)*(1-(b+1)*x),label=r'$\vec{v}_1$',linestyle='-.')
plt.plot(x,b/(a*x),label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = r*x+y+x^2$, $\dot{y} = -x+r*y+2*x^2$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.legend()

plt.show()
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]::
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
    
dt = 0.01
tf = 30

x, y = np.meshgrid(np.linspace(0, tf, 20), np.linspace(-2*np.pi, 2*np.pi, 30))
t = np.arange(0, tf, dt)

def drdt(r,theta,t): return 1

def dthetadt(r,theta,t): return np.sin(r) - np.sin(theta)

points = np.linspace(-2*np.pi, 2*np.pi, 30)

X = drdt(x,y,t)
Y = dthetadt(x,y,t)

plt.figure(dpi=300)
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
plt.quiver(x, y, X, Y, color='grey', width=0.003)
for point in points:
    t, theta = n_coupled_rk4([drdt, dthetadt], [0,point], dt, tf)
    plt.plot(t, theta, color='tab:blue')
for point in points:
    t, theta = n_coupled_rk4([drdt, dthetadt], [30,point], -dt, -tf)
    plt.plot(t, theta, color='tab:blue')
#for point in points:
#    t, theta = n_coupled_rk4([drdt, dthetadt], [10,point], -dt, -tf)
#    plt.plot(t, theta, color='tab:blue')
plt.title(r'$\dot{\theta}=\sin t-\sin\theta$')
plt.ylabel(r'$\theta$')
plt.xlabel(r'$t$')
plt.xlim(0, tf)
plt.ylim(-2*np.pi, 2*np.pi)
#plt.legend()
plt.show()


#%% In class

def f(x): return r*x*(1-x)

r = 3.828

x = np.arange(0, np.pi, 0.01)

xn = [0.15]
yn = [0]

for n in range(100):
    xn.append(xn[-1])
    yn.append(f(xn[-1]))
    xn.append(yn[-1])
    yn.append(yn[-1])

plt.figure(dpi=300)
plt.plot(x, f(x), label=r'$x_{n+1}=rx(1-x)$', color='tab:blue')
#plt.plot(x, r*x)
plt.plot(x, x, color='tab:orange')
plt.plot(xn, yn, linewidth=1.25, linestyle='--', color='tab:green')
plt.scatter(xn[1::2], yn[1::2], marker='.', color='tab:green')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel(r'$x_n$')
plt.ylabel(r'$x_{n+1}$')
plt.legend()
plt.show()


#%% Alice and Bob

dt = 0.01
tf = 10

t = np.arange(0, tf, dt)

def dudt(u,v,w,t): return sigma*(v-u)

def dvdt(u,v,w,t): return r*u - v - 10*u*w

def dwdt(u,v,w,t): return 5*u*(v - b*w)

r = 1
sigma = 1
b = 1

plt.figure(dpi=300)
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
u, v, w = n_coupled_rk4([dudt, dvdt, dwdt], [1,1,1], dt, tf)
plt.plot(t, u[:-1], label=r'$u$')
plt.plot(t, v[:-1], label=r'$v$')
plt.plot(t, w[:-1], label=r'$w$')
plt.xlabel(r'$t$')
plt.xlim(0, tf)
plt.ylim(-1, 1)
plt.legend()
plt.show()

u += 1*np.random.rand(len(u))

def durdt(ur,vr,wr,t): return sigma*(vr-ur)

def dvrdt(ur,vr,wr,t): return r*u[np.where(np.arange(0, tf, dt))==t] - vr - \
    10*u[np.where(np.arange(0, tf, dt))==t]*wr

def dwrdt(ur,vr,wr,t): return 5*u[np.where(np.arange(0, tf, dt))==t]*(vr - b*wr)

plt.figure(dpi=300)
plt.plot([-10,10], [0,0], color='black', linewidth=0.75, linestyle='--')
plt.plot([0,0], [-10,10], color='black', linewidth=0.75, linestyle='--')
ur, vr, wr = n_coupled_rk4([durdt, dvdt, dwdt], [1,1,1], dt, tf)
plt.plot(t, ur[:-1], label=r'$u_r$')
plt.plot(t, vr[:-1], label=r'$v_r$')
plt.plot(t, wr[:-1], label=r'$w_r$')
plt.xlabel(r'$t$')
plt.xlim(0, tf)
plt.ylim(-1, 1)
plt.legend()
plt.show()



# In[9]:
# In[9]:
# In[9]:
# In[9]::
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]::
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]::
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]: