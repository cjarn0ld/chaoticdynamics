#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:48:24 2023

@author: charlesarnold
"""


# In[9]:
# In[9]:
# In[9]:
# In[9]:
import matplotlib.pyplot as plt

# Define the baker's map function
def baker(x, y, a):
    if x < 0.5:
        return 2*x, a*y
    else: 
        return 2*x-1, a*y+0.5


# Define the number of iterations and the initial values
n = 1000
x0, y0 = 0.001, 0.001
a = 0.2

# Generate the sequence of points by iteratively applying the map
points = []
x, y = x0, y0
for i in range(n):
    points.append((x, y))
    x, y = baker(x, y, a)

# Plot the sequence of points
plt.scatter(*zip(*points), c=range(n), cmap='coolwarm')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
# In[9]:
# In[9]:
# In[9]:
import matplotlib.pyplot as plt
import numpy as np

# Define the baker's map function
def baker(x, y, a):
    if x < 0.5:
        xn = 2*x
        yn = a*y
    else:
        xn = 2*x - 1
        yn = a*y+0.5
    return xn, yn

# Define the number of iterations and the initial values
n = 1000
x0, y0 = 0.01, 0.01
a = 0.2

# Generate the sequence of points by iteratively applying the map
x, y = x0, y0
xs, ys = [], []
for i in range(n):
    if i > 0:
        xs.append(x)
        ys.append(y)
    x, y = baker(x, y, a)

# Generate the 2D return map by plotting (xn, yn) over (xn+1, yn+1)
xn = xs[:-1]
xn1 = xs[1:]
yn = ys[:-1]
yn1 = ys[1:]
plt.scatter(xn, yn, c=xn1, cmap='coolwarm')
plt.xlabel('x_n')
plt.ylabel('y_n')
plt.colorbar()
plt.show()

plt.scatter(xn1, yn1, c=xn, cmap='coolwarm')
plt.xlabel('x_{n+1}')
plt.ylabel('y_{n+1}')
plt.colorbar()
plt.show()
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]:
# In[9]: