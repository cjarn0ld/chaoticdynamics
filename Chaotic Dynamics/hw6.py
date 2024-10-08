#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 21:29:15 2023

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
dt = 0.01
tf = 5

x, y = np.meshgrid(np.linspace(0, 3, 15), np.linspace(0, 3, 15))
t = np.arange(0, tf, dt)

def dxdt(x,y,t): return x*(3-2*x-y)

def dydt(x,y,t): return y*(2-x-y)

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(0.1,0.1), (0.5,0.5), (0.3,1.5), (2,2)]

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
plt.plot(-2*x+3,y,label=r'$\vec{v}_1$',linestyle='-.')
plt.plot(-x+2,y,label=r'$\vec{v}_2$',linestyle='-.')
plt.title(r'$\dot{x} = x(3-2x-y)$, $\dot{y} = y(2-x-y)$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(0, 3)
plt.ylim(0, 3)
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

#question 5    
dt = 0.01
tf = 5

x, y = np.meshgrid(np.linspace(-5, 5, 15), np.linspace(-5, 5, 15))
t = np.arange(0, tf, dt)

def dxdt(x,y,t): return y

def dydt(x,y,t): return x*cos(y)

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(3/2,1), (3,-1), (-2.5,1), (1,-1)]

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
plt.title(r'$\dot{x} = y$, $\dot{y} = xcos(y)$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.legend(loc='lower right')

plt.show()


# In[9]:

# In[9]:
    
dt = 0.01
tf = 5

#question4

x, y = np.meshgrid(np.linspace(-5, 5, 15), np.linspace(-5, 5, 15))
t = np.arange(0, tf, dt)

def dxdt(x,y,t): return x*y

def dydt(x,y,t): return x*y-y

X = dxdt(x,y,t)
Y = dydt(x,y,t)

points = [(1,1), (3,-1), (-2.5,1), (1,-1)]

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
plt.title(f'$\dot{x} = x*y$, $\dot{y} = *x*y-y*{l}$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.legend(loc='lower right')

plt.show()

# In[9]:[9]:

# In[9]:

# In[9]:   
    """
    Runge-Kutta 4 for three coupled differential equations

    Parameters
    ----------
    dxdt : function
        the first differential equation to be integrated using RK4.
    dydt : function
        the second differential equation to be integrated using RK4.
    dzdt : function
        the third differential equation to be integrated using RK4.
    x_0 : float
        the initial value of the first dependent variable.
    y_0 : float
        the initial value of the second dependent variable.
    z_0 : float
        the initial value of the third dependent variable.
    dt : TYPE
        DESCRIPTION.
    tf : TYPE
        DESCRIPTION.

    dt : float
        the step in the independent variable.
    tf : float
        the final value of the independent variable to be integrated to.

    Returns
    -------
    three element tuple
        three numpy arrays, each corresponding to one of the three fully integrated
        dependent variables.

    """

# In[9]:

# In[9]:
def tripple_rk4(dxdt, dydt, dzdt, x_0, y_0, z_0, dt, tf):

    
    def f(xyz,t):
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        return np.array([dxdt(x,y,z,t), dydt(x,y,z,t), dzdt(x,y,z,t)], float)
    
    t_arr = np.arange(0, tf, dt)
    X = []
    Y = []
    Z = []
    
    X.append(x_0)
    Y.append(y_0)
    Z.append(z_0)
    
    xyz = np.array([x_0, y_0, z_0], float)
    
    for t in t_arr:
        k1 = f(xyz, t) * dt
        k2 = f(xyz + 0.5*k1, t + 0.5*dt) * dt
        k3 = f(xyz + 0.5*k2, t + 0.5*dt) * dt
        k4 = f(xyz + k3, t + dt) * dt
        xyz += (k1 + 2*k2 + 2*k3 + k4)/6
        X.append(xyz[0])
        Y.append(xyz[1])
        Z.append(xyz[2])
    
    return np.array(X), np.array(Y), np.array(Z)

# In[9]:

# In[9]:
# set up paramaters relating to time
dt = 0.01                   # time step
tf = 5                      # final time
t = np.arange(0, tf, dt)    # time array

r = np.array([1,0,-1])

for i in range(len(r)):

    def dxdt(x,y,z,t): return r[i]*x*z
    
    def dydt(x,y,z,t): return r[i]*y*z
    
    def dzdt(x,y,z,t): return -r[i]*x*z-r[i]*y*z

    x, y, z = tripple_rk4(dxdt, dydt, dzdt, 0.2, 0.3, 0.5, dt, tf)
    
    plt.figure()
    plt.plot(t, x[:-1], label=r'$Right$')
    plt.plot(t, y[:-1], label=r'$Left$')
    plt.plot(t, z[:-1], label=r'$Center$')
    plt.title(f'$r={r[i]}$')
    plt.ylabel(r'Political Population')
    plt.xlabel(r'$t$')
    plt.ylim(0,1)
    plt.legend()
    plt.show()


# In[9]:

# In[9]:

# In[9]:[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]: