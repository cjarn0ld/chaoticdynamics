#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 13:11:00 2023

@author: charlesarnold
"""

#GENERATES LORENZ ATTRACTOR
import numpy as np
from numpy import pi,sin,exp
from scipy.integrate import odeint,RK45,solve_ivp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from IPython.display import clear_output
import time
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema




# In[9]:
#number 1
def lorentz(xt, yt, zt):
    x=[];y=[];z=[];RKt=[];tt=0;dt=0.01;
    tf=2000;xt=xt;yt=yt;zt=zt;  # CHANGE INITIAL CONDITIONS, RUNTIME
    sig=10;b=8/3;
    r=100  # CHANGE R
    RKt.append(tt);x.append(xt);y.append(yt);z.append(zt);
    
    while(tt<tf) :
      x1=xt;y1=yt;z1=zt;
      xk1=sig*(y1-x1);yk1=r*x1-y1-x1*z1;zk1=x1*y1-b*z1
      x2=xt+dt*xk1/2;y2=yt+dt*yk1/2;z2=zt+dt*zk1/2;
      xk2=sig*(y2-x2);yk2=r*x2-y2-x2*z2;zk2=x2*y2-b*z2
      x3=xt+dt*xk2/2;y3=yt+dt*yk2/2;z3=zt+dt*zk2/2;
      xk3=sig*(y3-x3);yk3=r*x3-y3-x3*z3;zk3=x3*y3-b*z3
      x4=xt+dt*xk3;y4=yt+dt*yk3;z4=zt+dt*zk3;
      xk4=sig*(y4-x4);yk4=r*x4-y4-x4*z4;zk4=x4*y4-b*z4
      vx=(xk1+2*xk2+2*xk3+xk4)/6;vy=(yk1+2*yk2+2*yk3+yk4)/6;vz=(zk1+2*zk2+2*zk3+zk4)/6;
      xt=xt+dt*vx;yt=yt+dt*vy;zt=zt+dt*vz;tt=tt+dt
      x.append(xt);y.append(yt);z.append(zt);RKt.append(tt)
    
    
    fig=plt.figure(figsize=[10,10]);ax=fig.add_subplot(3,3,1);ax.set_xlabel('Time (t)');ax.set_ylabel('x(t)')
    ax.plot(RKt,x);ax=fig.add_subplot(3,3,2);ax.set_xlabel('Time (t)');ax.set_ylabel('y(t)')
    ax.plot(RKt,y);ax=fig.add_subplot(3,3,3);ax.set_xlabel('x(t)');ax.set_ylabel('z(t)')
    ax.plot(x,z);ax=fig.add_subplot(2,1,2,projection='3d');ax.set_xlabel('x');ax.set_ylabel('y');ax.set_zlabel('z');ax.plot3D(x,y,z)

# In[9]:
lorentz(1,1,1)
lorentz(0,5,1)
lorentz(1,0,5)

# In[9]:
#number 2

# Define the Lorenz system
def lorenz(t, state, sigma=10.0, beta=8/3, rho=28.0):
    x, y, z = state
    xdot = sigma * (y - x)
    ydot = x * (rho - z) - y
    zdot = x * y - beta * z
    return [xdot, ydot, zdot]

# Set the integration parameters
dt = 0.01
t0, tf = 0, 100
t = np.arange(t0, tf, dt)

# Set the initial conditions
state0 = [1.0, 1.0, 1.0]

# Integrate the Lorenz system
from scipy.integrate import solve_ivp
sol = solve_ivp(lorenz, [t0, tf], state0, t_eval=t)

# Compute the local maxima of z(t)
z_max = []
z_prev = sol.y[2, 0]
for z in sol.y[2, 1:]:
    if z_prev > z and len(z_max) < len(sol.t) - 1:
        z_max.append(z_prev)
    z_prev = z

# Compute the Lorenz map
z_max = np.array(z_max)
z_map = z_max[1:]
z_max = z_max[:-1]

# Plot the Lorenz map
plt.plot(z_max, z_map, 'o-', markersize=1)
plt.xlabel('$z_n$')
plt.ylabel('$z_{n+1}$')
plt.title('Lorenz Map')
plt.show()
print(f"The last value of z is {z[-1]}")
# In[9]:
    #number 3
# Define the Lorenz system
def lorenz(x, y, z, s=10, r=28, b=8/3):
    dxdt = s*(y - x)
    dydt = r*x - y - x*z
    dzdt = x*y - b*z
    return dxdt, dydt, dzdt

# Define the Lorenz map
def lorenz_map(xvals, yvals, zvals):
    x_max = []
    y_max = []
    z_max = []
    for i in range(len(xvals)-2):
        if yvals[i+1] > yvals[i] and yvals[i+1] > yvals[i+2]:
            x_max.append(xvals[i+1])
            y_max.append(yvals[i+1])
            z_max.append(zvals[i+1])
    xmap = np.array(x_max)
    ymap = np.array(y_max)
    zmap = np.array(z_max)
    return xmap, ymap, zmap

# Set initial conditions and time steps
x0, y0, z0 = (0.1 ,0.1, 0.1)
dt = 0.01
t = np.arange(0, 100, dt)

# Integrate the Lorenz system
x, y, z = np.zeros(len(t)), np.zeros(len(t)), np.zeros(len(t))
x[0], y[0], z[0] = x0, y0, z0
for i in range(len(t)-1):
    dxdt, dydt, dzdt = lorenz(x[i], y[i], z[i])
    x[i+1] = x[i] + dxdt * dt
    y[i+1] = y[i] + dydt * dt
    z[i+1] = z[i] + dzdt * dt

# Compute the Lorenz map
xmap, ymap, zmap = lorenz_map(x, y, z)

# Create cobweb plot
fig, ax = plt.subplots(figsize=(8,8))
ax.set_xlim(-20, 20)
ax.set_ylim(-20, 40)
ax.plot([-20, 20], [0, 0], 'k--')
ax.plot([0, 0], [-20, 40], 'k--')
ax.set_xlabel('z_n')
ax.set_ylabel('z_{n+1}')

# Plot the Lorenz map
ax.plot(xmap, ymap, markersize=3)

# Plot the cobweb
for i in range(len(xmap)-1):
    ax.plot([xmap[i], xmap[i]], [ymap[i], zmap[i]], 'r--', linewidth=0.5)
    ax.plot([xmap[i], xmap[i+1]], [zmap[i], ymap[i+1]], 'r-', linewidth=0.5)
ax.plot([xmap[-1], xmap[-1]], [ymap[-1], zmap[-1]], 'r--', linewidth=0.5)

plt.show()
# In[9]:

# Lorenz system parameters
sigma = 10
beta = 8/3
rho = 28

# Lorenz equations
def lorenz(x, y, z, sigma=sigma, beta=beta, rho=rho):
    x_dot = sigma * (y - x)
    y_dot = x * (rho - z) - y
    z_dot = x * y - beta * z
    return x_dot, y_dot, z_dot

# Parameters for numerical integration
dt = 0.01
num_steps = 10000

# Initial conditions
x, y, z = (1, 1, 1.05)

# Lists to store the Lorenz map
z_max = [z]
z_prev = z

# Numerical integration
for i in range(num_steps):
    # Compute derivatives
    x_dot, y_dot, z_dot = lorenz(x, y, z)

    # Update state
    x += x_dot * dt
    y += y_dot * dt
    z += z_dot * dt

    # Check for local maximum
    if z_dot < 0 and z_prev > z:
        z_max.append(z)
        z_prev = z

# Plot the Lorenz map
fig, ax = plt.subplots()
ax.plot(z_max[:-1], z_max[1:], 'bo', markersize=2)
ax.plot([0, 50], [0, 50], 'k--')
ax.set_xlabel('$z_n$')
ax.set_ylabel('$z_{n+1}$')
ax.set_title('Lorenz Map')

# Plot the cobweb map
fig2, ax2 = plt.subplots()
ax2.plot(z_max, z_max, 'bo', markersize=2)
x = 0
for i in range(len(z_max)-1):
    x = z_max[i]
    y = z_max[i+1]
    ax2.plot([x, y], [y, y], 'r')
    ax2.plot([y, y], [y, x], 'r')

ax2.plot([0, 50], [0, 50], 'k--')
ax2.set_xlabel('$z_n$')
ax2.set_ylabel('$z_{n+1}$')
ax2.set_title('Cobweb Map')
plt.show()

# In[9]:# In[9]:
#number 4 

def lorentz1(r):
    xt = 0.1
    yt = 0.1 
    zt = 0.1
    x=[];y=[];z=[];RKt=[];tt=0;dt=0.01;
    tf=2000;xt=xt;yt=yt;zt=zt;  # CHANGE INITIAL CONDITIONS, RUNTIME
    sig=10;b=8/3;
    r=r # CHANGE R
    RKt.append(tt);x.append(xt);y.append(yt);z.append(zt);
    
    while(tt<tf) :
      x1=xt;y1=yt;z1=zt;
      xk1=sig*(y1-x1);yk1=r*x1-y1-x1*z1;zk1=x1*y1-b*z1
      x2=xt+dt*xk1/2;y2=yt+dt*yk1/2;z2=zt+dt*zk1/2;
      xk2=sig*(y2-x2);yk2=r*x2-y2-x2*z2;zk2=x2*y2-b*z2
      x3=xt+dt*xk2/2;y3=yt+dt*yk2/2;z3=zt+dt*zk2/2;
      xk3=sig*(y3-x3);yk3=r*x3-y3-x3*z3;zk3=x3*y3-b*z3
      x4=xt+dt*xk3;y4=yt+dt*yk3;z4=zt+dt*zk3;
      xk4=sig*(y4-x4);yk4=r*x4-y4-x4*z4;zk4=x4*y4-b*z4
      vx=(xk1+2*xk2+2*xk3+xk4)/6;vy=(yk1+2*yk2+2*yk3+yk4)/6;vz=(zk1+2*zk2+2*zk3+zk4)/6;
      xt=xt+dt*vx;yt=yt+dt*vy;zt=zt+dt*vz;tt=tt+dt
      x.append(xt);y.append(yt);z.append(zt);RKt.append(tt)
    
    
    fig=plt.figure(figsize=[10,10]);ax=fig.add_subplot(3,3,1);ax.set_xlabel('Time (t)');ax.set_ylabel('x(t)')
    ax.plot(RKt,x);ax=fig.add_subplot(3,3,2);ax.set_xlabel('Time (t)');ax.set_ylabel('y(t)')
    ax.plot(RKt,y);ax=fig.add_subplot(3,3,3);ax.set_xlabel('x(t)');ax.set_ylabel('z(t)')
    ax.plot(x,z);ax=fig.add_subplot(2,1,2,projection='3d');ax.set_xlabel('x');ax.set_ylabel('y');ax.set_zlabel('z');ax.plot3D(x,y,z)

# In[9]:
lorentz1(166.3)
lorentz1(212)
lorentz1(156)
# In[9]:
    
def xyz_rk4(dxdt, dydt, dzdt, x0, y0, z0, dt, tf):
    
    def xyz(V,t):
        x = V[0]
        y = V[1]
        z = V[2]
        xyz = np.array([dxdt(x,y,z,t), dydt(x,y,z,t), dzdt(x,y,z,t)], float)
        return xyz
    
    T = np.arange(0, tf, dt)
    X = []
    Y = []
    Z = []
    V = []
    
    X.append(x0)
    Y.append(y0)
    Z.append(z0)
    
    V = np.array([x0, y0, z0], float)
    
    for i in T:
        k1 = xyz(V, t) * dt
        k2 = xyz(V + 0.5*k1, t + 0.5*dt) * dt
        k3 = xyz(V + 0.5*k2, t + 0.5*dt) * dt
        k4 = xyz(V + k3, t + dt) * dt
        V = V+(k1 + 2*k2 + 2*k3 + k4)/6
        X.append(V[0])
        Y.append(V[1])
        Z.append(V[2])
    
     
    x = np.array(X)
    y = np.array(Y) 
    z = np.array(Z)
    return x,y,z

# In[9]:

    
dt = 0.001
tf = 100

t = np.arange(0, tf, dt)

def x_dot(x,y,z,t): 
    return sigma * (y - x)

def y_dot(x,y,z,t): 
    return r * x - y - x * z

def z_dot(x,y,z,t): 
    return x * y - b * z

sigma = 10
r = 28
b = 8/3

x, y, z = xyz_rk4(x_dot, y_dot, z_dot,1,1,1, dt, tf)

zn = z[argrelextrema(z, np.greater)[0]]

xt = np.linspace(np.min(zn),np.max(zn),len(zn))
zn, zn1 = zn[0:-2], zn[1:-1]
yt = xt

def cobweb(x1, y1, znarray, zn1array, N):
    #setup
    x = x1 # xi
    y = y1 # yi
    x_arr = []
    y_arr = []
    x_arr.append(x)
    y_arr.append(y) # append starting point

    # cobweb
    for i in range(N):
        # find y value from conditions
        y_prime = np.max(np.extract(np.abs(zn - float(x)) < 0.1, zn1))
        # append new points
        y_arr.append(float(y_prime))
        x_arr.append(float(x))

        # set y equal to x
        x = y_prime
        y_arr.append(float(y_prime))
        x_arr.append(float(x))

    return x_arr, y_arr

x, y = cobweb(37, 20, zn, zn1, 1000)

plt.figure()
plt.scatter(zn,zn1, color='tab:blue')
plt.plot(x, y, color='tab:purple')
plt.plot([29,48], [29,48], color='tab:orange')
plt.scatter(x[1::2], y[1::2], marker='.', color='black')
plt.ylim(27,46)

# In[9]:
# In[9]:
# In[9]:

# In[9]:# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]: