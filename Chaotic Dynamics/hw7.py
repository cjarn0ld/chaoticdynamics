#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:20:41 2023

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
dt = 0.01
tf = 10

x, y = np.meshgrid(np.linspace(-3, 3, 15), np.linspace(-3, 3, 15))
t = np.arange(0, tf, dt)

a = 1

def dxdt(x,y,t): return y

def dydt(x,y,t): return -a*y*(x**2+y**2-1)-x

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(1,2), (0.2,-0.2), (-2,1), (-0.5,-0.25)]

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
#plt.plot(-2*x+3,y,label=r'$\vec{v}_1$',linestyle='-.')
#plt.plot(-x+2,y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = y$, $\dot{y} =  -ay(x^2+y^2-1)-x$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.legend(loc='lower right')

plt.show()


# In[9]:
# In[9]:dt = 0.01
tf = 10

x, y = np.meshgrid(np.linspace(-5, 5, 15), np.linspace(-5, 5, 15))
t = np.arange(0, tf, dt)

a = 1

def dxdt(x,y,t): return y+2*x*y

def dydt(x,y,t): return x+x**2-y**2

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(0.1,0.1), (0,4), (-3,1), (-4,-3)]

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
#plt.plot(-2*x+3,y,label=r'$\vec{v}_1$',linestyle='-.')
#plt.plot(-x+2,y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = y+2xy$, $\dot{y} =  -x+x^2-y^2x$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.legend(loc='lower right')

plt.show()
# In[9]:
# In[9]:
tf = 10

x, y = np.meshgrid(np.linspace(-5, 5, 15), np.linspace(-5, 5, 15))
t = np.arange(0, tf, dt)

a = 1

def dxdt(x,y,t): return x-y-x*(x**2+5*y**2) #((x**2+y**2)-(x**4+y**4)-6*x**2*y**2)/(x**2+y**2)**(0.5)

def dydt(x,y,t): return x+y-y*(x**2+y**2) #(x**2+y**2+4*x*y**3)/(x**2+y**2)

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(0.1,0.1), (0,4), (-3,1), (-4,-3)]

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
#plt.plot(-2*x+3,y,label=r'$\vec{v}_1$',linestyle='-.')
#plt.plot(-x+2,y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = y+2xy$, $\dot{y} =  -x+x^2-y^2x$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.legend(loc='lower right')

plt.show()
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
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:v