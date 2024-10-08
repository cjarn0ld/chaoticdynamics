#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 12:36:05 2023

@author: charlesarnold
"""

print('Homework #1, due Sunday Jan 20')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import math
from scipy.signal import find_peaks
from scipy import signal
from PIL import Image


# In[9]:


# In[9]:
print("2.7.1 Potentials")
print("For each of the following vector fields, plot the potential function V(x) and identify all the equilibrium points and their stability.")
print(r'$ \dot{x} $ =(1−x)')
# In[9]:

    
# In[9]:
print("3.1 Saddle-Node Bifurcation") 
print("For each of the following exercises, sketch all the qualitatively different vector fields that occur as r is varied. Show that a saddle-node bifurcation occurs at a critical value of r, to be determined. Finally, sketch the bifurcation diagram of fixed points x* versus r.") 

print("3.1.2") 
print("xdot = r−cosh((x)") 

print("3.2.4") 
print("xdot=x(r−e^x)")

# In[9]:

# In[9]:
print("3.4 Pitchfork Bifurcation")
print("In the following exercises, sketch all the qualitatively different vector fields that occur as r is varied. Show that a pitchfork bifurcation occurs at a critical value of r (to be determined) and classify the bifurcation as supercritical or subcritical. Finally, sketch the bifurcation diagram of x* vs. r.")
print("The next exercises are designed to test your ability to distinguish among the various types of bifurcations—it’s easy to confuse them! In each case, find the values of r at which bifurcations occur, and classify those as saddle-node, transcritical, supercritical pitchfork, or subcritical pitchfork. Finally, sketch the bifurcation dia- gram of fixed points x* vs. r.")

print("3.4.6")
print("x = rx−x/1+x ")
# In[9]:
print("3.4.14")
print("(Subcritical pitchfork) Consider the system x = rx + x3 − x5 , which exhibits a subcritical pitchfork bifurcation. ")
print("a) Find algebraic expressions for all the fixed points as r varies")
print("")
print("b)  Sketch the vector fields as r varies. Be sure to indicate all the fixed points and their stability.")
print("")
print("c)  Calculate rs , the parameter value at which the nonzero fixed points are born in a saddle-node bifurcation")
# In[9]:
print("	6. (a) Write a computer program to integrate the logistic equation ̇N = rN (1 − N/k) for several values of r (both positive and negative) and k and put all the curves on a single plot.")
print("")
print("b) Now, write a “wrapper” program to identify the stable fixed point as a function of r like this. Beginning with r = −1 and N != 0")
print("i. run your program from part (a) to a final time long enough for the system to have converged to the fixed point")
print("ii. Record the value of r and the last value from your program")
print("iii. increase r by a small amount (0.01? 0.1?)")
print(" Repeat until r = 1 and then plot your “last value” on the y-axis and the corre-sponding value of r on the x − axis to show the transfer of stability!")
print("")
print("(c) Bonus: Explore how your plot from (b) changes as you vary tf , the final time to which you integrate in step (b)(i). In particular, what happens when r is near 0? How can we explain this behavior?")



# In[9]:

 


def integral(N,r,k):
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

r = np.linspace(-1,3,5)
k = np.linspace(3,5,5)

c = np.linspace(1,5,5)
for j in range(len(c)):
    plt.figure()
    for i in range(len(k)):
        plt.plot(*integral(round(c[j]),r[i],k[i]), linestyle = "-", label = f'N = {round(c[j],3)}; r = {r[i]}; k = {k[i]}') 
        plt.xlabel("Time") 
        plt.ylabel("N") 
        plt.legend()
        plt.title("my code wants to break :( ")  
    




# In[9]:
r1 = np.linspace(-1,5,1000)
y = []
for i in range(len(r1)):
    y.append(integral(2,r1[i],3)[1][-1])
plt.plot(r1, y)
plt.title('r vs xf')
plt.xlabel('r')
plt.ylabel('xf')
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

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:
# In[9]:
