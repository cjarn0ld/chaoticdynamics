#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 14:59:52 2023

@author: charlesarnold
"""
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
import numpy as np
from numpy import sqrt
import random
from scipy.integrate import odeint,RK45,solve_ivp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from numpy import sqrt
import random
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import log,exp
from scipy import stats


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

def gen_points(start, end, n_points):
    points = []
    for i in np.linspace(start,end,n_points):
        points.append((i,0))
        points.append((0,i))
        points.append((i,i))
        points.append((i,-i))
    return points
# In[9]:
dt = 0.01
tf = 50

#question4

x, y = np.meshgrid(np.linspace(-5, 5, 14), np.linspace(-5, 5, 14))
t = np.arange(0, tf, dt)

a = 2
b = 2

def dxdt(x,y,t): return a-x+x**2*y

def dydt(x,y,t): return b-x**2*y

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
    line = plt.plot(*coupled_rk4(dxdt, dydt, point[0], point[1], dt, tf), color='tab:blue')
    plt.plot(*coupled_rk4(dxdt, dydt, point[0], point[1], -dt, -tf), color=line[0].get_color())
    
x = np.arange(-10,10,0.1)
y = x

plt.plot(x,(x-a)/x**2,label=r'$\dot{x} = 0$', color = 'purple', linestyle='-.')
plt.plot(x,(b/x**2),label=r'$\dot{y} = 0$',linestyle='-.')
plt.title(r'$\dot{x} = r*x+y+x^2$, $\dot{y} = -x+r*y+2*x^2$ $, ' f' a = {a}, b = {b}$')
plt.ylabel(r'$y$')
plt.xlabel(r'$x$')
plt.xlim(0, 5)
plt.ylim(-1, 5)



plt.show()

# In[9]:
# In[9]:
# In[9]:

# In[9]:
    #number 2
R_list = np.linspace(0, 4, 1000)
x0 = 0.3
N = 500 #number of iterations

def logis(r):
    x_list = [x0]  #make x0  
    for i in range(N-1):
        x_list.append(r * x_list[-1] - (x_list[-1])**3) #itterate cubic function
    return x_list[400:]

#make the empty arrays 
x = []
R = []
for r in R_list:
    x.append(logis(r))
    R.append([r] * 100) 
    
x = np.array(x).ravel() #flatten the array for x and r for each, get rid of the extra points for each x
R = np.array(R).ravel() 


plt.figure(figsize=(11, 8))
plt.xlabel('r',fontsize=20)
plt.ylabel(r'$x_n^{*}$',fontsize=20)
plt.title(f'\nThe bifurcation diagram of the Cubic Map\n\n 0< R < 3 | 'r'$x_{n+1}=rx_{n}-x_n^3$|  'r'$x_0$'f'={x0}',fontsize=20)
plt.scatter(R, x,marker='.', color='tab:blue', s=0.1)
plt.xlim(0.8,3)
plt.show()


# In[9]:
    #number 3
a = 1.4
b =0.3
N = 100000# number of iterations

x = np.zeros(N)
y = np.zeros(N)


# initial condition
x[0] = 0.4
y[0] = 0.3

for i in range(N-1):
    x[i+1] = 1 + y[i] - a * abs(x[i])
    y[i+1] = b * x[i]
        
# plotting the attractor    
plt.figure(figsize=(11, 8))
plt.scatter(x, y)#, '.', markersize=0.5)
plt.title(f'Strange attractor for "Henon like Map" with a={a},b={b}',fontsize = 20)
plt.xlabel('$x_n$',fontsize = 20)
plt.ylabel('$y_n$',fontsize = 20)
plt.show()


# In[9]:
a = 1.4
b = np.linspace(-1,2,4)
N = 1000 # number of iterations

x = np.zeros(N)
y = np.zeros(N)


# initial condition
x[0] = 0.4
y[0] = 0.3

for i in range(len(b)):
    # iterating the map    

    for j in range(N-1):
        x[j+1] = 1 + y[j] - a * abs(x[j])
        y[j+1] = b[i] * x[j]
        
    # plotting the attractor    
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.set_title(f'Strange attractor for for "Henon Map" with a={a},b={b[i]}')
    ax.set_xlabel('$x_n$',fontsize = 20)
    ax.set_ylabel('$y_n$',fontsize = 20)
    ax.plot(x, y)



# In[9]:

#make histagram bins to find box dim and coralationdim

n=len(x)
print(n)
d=[];
i=random.randint(0,n-1)
j=0
while (j<n):
  d.append(sqrt((x[i]-x[j])**2+(y[i]-y[j])**2))
  j+=1
d=sorted(d)
#plt.loglog(d,range(n))

xmin=min(min(x),min(y))
xmax=max(max(x),max(y))
L=xmax-xmin
edges=np.linspace(min(min(x),min(y)),max(max(x),max(y)),num=2**13+1)

(h,xedge,yedge,image)=plt.hist2d(x,y,bins=edges)
print(h.shape)


# In[9]:

e=[];N=[];
n=int(log(len(h))/log(2))-1
h2=h;
e.append(len(h2)/L)
N.append(np.count_nonzero(h2))
while(n>=0):
  h2=h2.reshape((2**n,2,2**n,-1)).sum(axis=3).sum(1)
  e.append(len(h2)/L)
  N.append(np.count_nonzero(h2))
  n-=1

pnty=range(1,len(d))  
box=np.polyfit(log(e),log(N),1)
pntdim=np.polyfit(log(d[1:]),log(pnty),1)

box_e = [exp(box[1])*i**box[0] for i in e]
pnt_e= [exp(pntdim[1])*i**pntdim[0] for i in d]

print("Box dimension=",box[0]," Pointwise dimension=",pntdim[0])
slope, intercept, r_value, p_value, std_err = stats.linregress(e, box_e)

plt.yscale('log')
plt.xscale('log')
plt.scatter(e,N)
plt.plot(e,box_e)
plt.title("Correlation Dimension Fit")
print("Slope: ", slope)
print("Intercept: ", intercept)
print("R-value: ", r_value)
print("P-value: ", p_value)
print("Standard error: ", std_err)
# In[9]:
#question 4

x=[];y=[];z=[];RKt=[];tt=5;dt=0.01;
tf=20;xt=2;yt=2;zt=2;  # CHANGE INITIAL CONDITIONS, RUNTIME
a=0.4;b=2;
c=4  # CHANGE R
RKt.append(tt);x.append(xt);y.append(yt);z.append(zt);
fig=plt.figure(figsize=[10,10]);ax=fig.add_subplot(3,3,1);ax.set_xlabel('Time (t)');ax.set_ylabel('x(t)')
plt.title(f"IC:(x,y,z) = ({xt},{yt},{zt}) t0 = {tt}")   
while(tt<tf) :
     x1=xt;y1=yt;z1=zt;
     xk1=(-y1-z1);yk1=(x1+a*y1);zk1=b+z1*(x1-c)
     x2=xt+dt*xk1/2;y2=yt+dt*yk1/2;z2=zt+dt*zk1/2;
     xk2=(-y2-z2);yk2=x2+a*y2;zk2=b+z2*(x2-c)
     x3=xt+dt*xk2/2;y3=yt+dt*yk2/2;z3=zt+dt*zk2/2;
     xk3=(-y3-z3);yk3=x3+a*y3;zk3=b+z3*(x3-c)
     x4=xt+dt*xk3;y4=yt+dt*yk3;z4=zt+dt*zk3;
     xk4=(-y4-z4);yk4=x4+a*y4;zk4=b+z4*(x4-c)
     vx=(xk1+2*xk2+2*xk3+xk4)/6;vy=(yk1+2*yk2+2*yk3+yk4)/6;vz=(zk1+2*zk2+2*zk3+zk4)/6;
     xt=xt+dt*vx;yt=yt+dt*vy;zt=zt+dt*vz;tt=tt+dt
     x.append(xt);y.append(yt);z.append(zt);RKt.append(tt)
   
     
    
    
ax.plot(RKt,x);ax=fig.add_subplot(3,3,2);ax.set_xlabel('Time (t)');ax.set_ylabel('y(t)')
ax.plot(RKt,y);ax=fig.add_subplot(3,3,3);ax.set_xlabel('Time (t)');ax.set_ylabel('z(t)')
ax.plot(RKt,z);ax=fig.add_subplot(2,1,2,projection='3d');ax.set_xlabel('x');ax.set_ylabel('y');ax.set_zlabel('z');ax.plot3D(x,y,z)

# In[9]:
lorentz(2,2,2)
lorentz(0,-5,1)
lorentz(-1,0,5)


# In[9]:

n=len(x)
print(n)
d=[];
i=random.randint(0,n-1)
j=0
while (j<n):
  d.append(sqrt((x[i]-x[j])**2+(y[i]-y[j])**2))
  j+=1
d=sorted(d)
#plt.loglog(d,range(n))

xmin=min(min(x),min(y))
xmax=max(max(x),max(y))
L=xmax-xmin
edges=np.linspace(min(min(x),min(y)),max(max(x),max(y)),num=2**13+1)

(h,xedge,yedge,image)=plt.hist2d(x,y,bins=edges)
print(h.shape)


# In[9]:

e=[];N=[];
n=int(log(len(h))/log(2))-1
h2=h;
e.append(len(h2)/L)
N.append(np.count_nonzero(h2))
while(n>=0):
  h2=h2.reshape((2**n,2,2**n,-1)).sum(axis=3).sum(1)
  e.append(len(h2)/L)
  N.append(np.count_nonzero(h2))
  n-=1

pnty=range(1,len(d))  
box=np.polyfit(log(e),log(N),1)
pntdim=np.polyfit(log(d[1:]),log(pnty),1)

box_e = [exp(box[1])*i**box[0] for i in e]
pnt_e= [exp(pntdim[1])*i**pntdim[0] for i in d]

print("Box dimension=",box[0]," Pointwise dimension=",pntdim[0])
coefficients = np.polyfit(e, N, 10)
# Evaluate the fitted polynomial

y_fit = np.polyval(coefficients,e)
# Print the coefficients
print("Coefficients:", coefficients)


print(slope)
plt.scatter(e, N)
plt.plot(e, box_e, label="Box Dimension")
plt.plot(e, y_fit, label='Polynomial Fit of degree 10')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("log(epsilon)")
plt.ylabel("log(N)")
plt.title("Correlation Dimension Fit")
plt.legend()
plt.show()
# In[9]:
def lorentz2(xt, yt, zt,tau):
    x=[];y=[];z=[];RKt1=[];RKt2=[];x2=[];tt=5;dt=0.01;
    tf=20;xt0=xt;yt0=yt;zt0=zt;  # CHANGE INITIAL CONDITIONS, RUNTIME
    a=0.4;b=2;
    c=4  # CHANGE R
    RKt1.append(tt);x.append(xt0);y.append(yt0);z.append(zt0)

    
    while(tt<tf) :
      x1=xt;y1=yt;z1=zt;
      xk1=(-y1-z1);yk1=(x1+a*y1);zk1=b+z1*(x1-c)
      x2=xt+dt*xk1/2;y2=yt+dt*yk1/2;z2=zt+dt*zk1/2;
      xk2=(-y2-z2);yk2=x2+a*y2;zk2=b+z2*(x2-c)
      x3=xt+dt*xk2/2;y3=yt+dt*yk2/2;z3=zt+dt*zk2/2;
      xk3=(-y3-z3);yk3=x3+a*y3;zk3=b+z3*(x3-c)
      x4=xt+dt*xk3;y4=yt+dt*yk3;z4=zt+dt*zk3;
      xk4=(-y4-z4);yk4=x4+a*y4;zk4=b+z4*(x4-c)
      vx=(xk1+2*xk2+2*xk3+xk4)/6;vy=(yk1+2*yk2+2*yk3+yk4)/6;vz=(zk1+2*zk2+2*zk3+zk4)/6;
      xt1=xt+dt*vx;yt=yt+dt*vy;zt=zt+dt*vz;tt1=tt+dt;tt2=tt+dt
      x.append(xt1);y.append(yt);z.append(zt);RKt1.append(tt1)
    
    fig=plt.figure(figsize=[10,10])
    ax.set_xlabel('Time (t)')
    ax.set_ylabel('x(t)')
    plt.title(f'IC:(x,y,z) = ({xt},{yt},{zt}, 'r'$\uptau$'f'= {tau}')
    ax.plot(RKt1,x)
    
# In[9]:
lorentz2(2,2,2,2)
lorentz2(0,-5,1,2)
lorentz2(-1,0,5,2)
# In[9]:
def lorentz2(xt, yt, zt, tau):
    x = [xt]; y = [yt]; z = [zt]; t = [0]
    a = 0.4; b = 2; c = 4; dt = 0.01; tf = 20
    
    while t[-1] < tf:
        xt1 = x[-1]; yt1 = y[-1]; zt1 = z[-1]
        xk1 = -yt1 - zt1
        yk1 = xt1 + a * yt1
        zk1 = b + zt1 * (xt1 - c)
        xt2 = xt1 + dt * xk1 / 2
        yt2 = yt1 + dt * yk1 / 2
        zt2 = zt1 + dt * zk1 / 2
        xk2 = -yt2 - zt2
        yk2 = xt2 + a * yt2
        zk2 = b + zt2 * (xt2 - c)
        xt3 = xt1 + dt * xk2 / 2
        yt3 = yt1 + dt * yk2 / 2
        zt3 = zt1 + dt * zk2 / 2
        xk3 = -yt3 - zt3
        yk3 = xt3 + a * yt3
        zk3 = b + zt3 * (xt3 - c)
        xt4 = xt1 + dt * xk3
        yt4 = yt1 + dt * yk3
        zt4 = zt1 + dt * zk3
        xk4 = -yt4 - zt4
        yk4 = xt4 + a * yt4
        zk4 = b + zt4 * (xt4 - c)
        x.append(xt1 + dt * (xk1 + 2 * xk2 + 2 * xk3 + xk4) / 6)
        y.append(yt1 + dt * (yk1 + 2 * yk2 + 2 * yk3 + yk4) / 6)
        z.append(zt1 + dt * (zk1 + 2 * zk2 + 2 * zk3 + zk4) / 6)
        t.append(t[-1] + dt)
    
    # Plot reconstructed attractor
    fig, ax = plt.subplots()
    ax.scatter(x[:-tau], x[tau:], s=1)
    ax.set_xlabel('x(t)')
    ax.set_ylabel('x(t + tau)')
    ax.set_title('Reconstructed attractor (tau = {})'.format(tau))
    
lorentz2(2,2,2,2)  
# In[9]:

n=len(x)
print(n)
d=[];
i=random.randint(0,n-1)
j=0
while (j<n):
  d.append(sqrt((x[i]-x[j])**2+(y[i]-y[j])**2))
  j+=1
d=sorted(d)
#plt.loglog(d,range(n))

xmin=min(min(x),min(y))
xmax=max(max(x),max(y))
L=xmax-xmin
edges=np.linspace(min(min(x),min(y)),max(max(x),max(y)),num=2**13+1)

(h,xedge,yedge,image)=plt.hist2d(x,y,bins=edges)
print(h.shape)

  
# In[9]:

# In[9]:
# In[9]:
  
e=[];N=[];
n=int(log(len(h))/log(2))-1
h2=h;
e.append(len(h2)/L)
N.append(np.count_nonzero(h2))
while(n>=0):
  h2=h2.reshape((2**n,2,2**n,-1)).sum(axis=3).sum(1)
  e.append(len(h2)/L)
  N.append(np.count_nonzero(h2))
  n-=1

pnty=range(1,len(d))  
box=np.polyfit(log(e),log(N),1)
pntdim=np.polyfit(log(d[1:]),log(pnty),1)

box_e = [exp(box[1])*i**box[0] for i in e]
pnt_e= [exp(pntdim[1])*i**pntdim[0] for i in d]

print("Box dimension=",box[0]," Pointwise dimension=",pntdim[0])
slope, intercept, r_value, p_value, std_err = stats.linregress(e, box_e)

plt.yscale('log')
plt.xscale('log')
plt.scatter(e,N)
plt.plot(e,box_e)
plt.title("Correlation Dimension Fit")
print("Slope: ", slope)
print("Intercept: ", intercept)
print("R-value: ", r_value)
print("P-value: ", p_value)
print("Standard error: ", std_err)

# In[9]:
# Calculate box counting dimension

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

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]:

# In[9]: