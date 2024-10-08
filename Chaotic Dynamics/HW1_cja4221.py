#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:27:41 2023

@author: charlesarnold
"""
print('Homework #1, due Friday Jan 20')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy.signal import find_peaks
from scipy import signal
from PIL import Image



# In[9]:

print('1. 2.2.4')
def equ1(x):
    x_dot = np.exp(-x)*np.sin(x)
    return x_dot
y=np.zeros(1000)
x = np.linspace(-10,10,1000)
plt.plot(x,equ1(x))
plt.plot(x,y)
plt.title(f"e^-x*sin(x)")
plt.xlabel("X")
plt.ylabel("X_dot")
plt.xlim(-10,-1)
plt.ylim(-200,200)

print("When finding stability points, setting sin(theta) = 0 gives us values of n*pi. this means that values of n*pi are local stable or unstable.")
print("when looking at the phase plot, when n*pi is odd it is a local stable point. whe n*pi is even it becomes a local unstable point.")


# In[9]:
print('2. 2.3.3')
print("2.3.3 ")
print("(Tumor growth) The growth of cancerous tumors can be modeled by the Gompertz law")
print("Ndot= -aN*ln(bN), where N ( t ) is proportional to the number of cells in the tumor, ")
print("and a, b > 0 are parameters.")




print("a) Interpret a and b biologically.")
print("a is the constant that determines the death rate of a cancer population while b is the growth rate of the cell population")
print("")
print("b) Sketch the vector field and then graph N ( t ) for various initial values.")
print("The predictions of this simple model agree surprisingly well with data on tumor")
print("growth, as long as N is not too small; see Aroesty et al. (1973) and Newton (1980) for examples.")
def equ3(N):
    n_dot = -N*np.log(N)
    return n_dot
dot = '\u0307'
N = np.linspace(0,10,1000)
plt.figure()
plt.plot(N, equ3(N))
plt.xlabel('b$\cdot$N')
plt.ylabel(r'$\dfrac{b}{a}$ $\cdot$ $ \dot{N} $')
plt.xlim(0,2)
plt.ylim(-2,4)
plt.plot(N, np.zeros(1000))
plt.title("-aN$\cdot$ln(b$\cdot$N)")

#desmose slider for a and b
#do linear stability with new fumction 

# In[9]:
print('3. 2.4.5')
print("")
print("Use linear stability analysis to classify the fixed points of the following systems. If linear stability analysis fails because f x′ ∗( ) , 0 use a graphical argument to decide the stability.")
print("")
def equ2(x):
    x_dot = 1-np.exp(-x**2)
    return x_dot
print("when plotting and doing linear stability, the derivative of x_dot shows that at the critical point the function is neither increasing or decreasing. Meaning that the function has a 0 slope at c.p. x=0.")
print("this means that to identify the stable points we have to graph it")
y=np.zeros(100)
x = np.linspace(-5,5,100)
plt.plot(x,equ2(x))
plt.plot(x,y)
plt.title(f"1-e^(-x^2)")
plt.xlabel("X")
plt.ylabel("X_dot")
print("after plotting, we find that the critical point(c.p.) @ x=0 is a semi stable point. when coming from the left side the flow is always decreasing so it will always approach 0. any position on the right, the graph is increaing and thus is an unstble point when an external force acts on the right, therefore making it a semi stable point")


#do linear stability in pucture 

# In[9]:
print(f"#{4} For the system N = N (r − a(N − b)^2)")
print(f"(a) Write a computer program to integrate the system (Euler is fine) for a particular a,b (your choice) and plot the behavior for a variety of initial values of N .")



def integral(N):
    a = 1    
    b = 2    
    r = 3
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
       N_dot = N_n*(r-a*(N_n-b)**2)
       time[i] = t + dt
       N_0[i] = N_n+N_dot*dt
      
    return time, N_0
c = np.linspace(1,8,7)
for i in range(1, len(c)+1):
    plt.plot(*integral(i), linestyle = "-", label = f'N = {i}') 
    plt.xlabel("Time")
    plt.ylabel("N")
    plt.title("a=1, b=2,r=3, dt=0.01")
    plt.ylim(0,7)
    plt.xlim(10,11)
    plt.legend()





# In[9]:
print(f"(b) Bonus: Find (analytically) the fixed points and their stability (for variable a, b)")
print("")
print(f"Plot (from your program) the approach toward stable fixed points (or away from unstable fixed points) on semi-log scales and graphically identify the Lyapunovexponent, comparing with your analytic solutions.")
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

# In[9]:

# In[9]:

# In[9]:
