# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:59:57 2019

@author: cami_
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import os
import time

sys.path.append(os.path.abspath("E:/JAVERIANA/COMPUTACION/code_proyecto2/RMS.py"))
from RMS import*

def dxdy(x,y,r):
   # f = float(-(2*x)*(y*y))
   f = r*y
   return(f)
    
def RK2(function,t0,y0,xmax,h):
    
    intervals = int(xmax/h)
    print("intervals:"+str(intervals))
    x = np.zeros((intervals,2))
    t=t0
    y=y0
    print("t:"+str(t)+ " | y:"+str(y))
    x[0,:] = [t,y]

    for i in range(1,intervals,1):
        
        k1 = h * dxdy(t,y,r)
        k2 = h * dxdy(t+h,y+k1,r)
        t = t + h
        y = y + 1/2*(k1+k2)
        print("t:"+str(t)+ " | y:"+str(y))
        print("k1:"+str(k1)+ " | k2:"+str(k2))
        x[i,:] = [t,y]
    return(x)
###################################################################
def growth(t0,N0,xmax,h,r):
    #tamano = muestras + 1
    intervals = int(xmax/h)
    print("intervals:"+str(intervals))
    malthus = np.zeros(shape=(intervals,2),dtype=float)
    malthus[0,0] = t0
    malthus[0,1] = N0
    
    for i in range(0,(intervals-1),1):
        malthus[i+1,0] = malthus[i,0] + h
        malthus[i+1,1] = N0*(math.e**(r*malthus[i+1,0]))
        #a*malthus[i,1]*(1-malthus[i,1]/K)*
        print(malthus[i+1,0])
    return(malthus);
    
    """   
        plt.plot(malthus[:,0],malthus[:,1],"ko") 
        plt.title("Malthus Growth")
        plt.xlabel("Seconds")
        plt.ylabel("Population")
        plt.show()
 """
#run program
##################################################################        
    
    
r=0.1   
X = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
analitics = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)
r=0.5
X2 = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
analitics2 = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)
r=1
X3 = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
analitics3 = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)
r=5
X4 = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
analitics4 = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)


"""
plt.plot(X.T[0],X.T[1],label='r=0.1 EST')
plt.plot(X2.T[0],X2.T[1],label='r=0.5 EST')
plt.plot(X3.T[0],X3.T[1],label='r=1 EST')
#plt.plot(X4.T[0],X4.T[1],label='r=5 EST')
"""
time_start = time.perf_counter()
r=0.1   
ff_list1 = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
time_elapsed1 = (time.perf_counter() - time_start)

time_start = time.perf_counter()
r=0.5   
ff_list2 = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
time_elapsed2 = (time.perf_counter() - time_start)

time_start = time.perf_counter()
r=1  
ff_list3 = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
time_elapsed3 = (time.perf_counter() - time_start)

r=5   
ff_list4 = RK2(function=dxdy,t0=0,y0=1E6,xmax=20,h=0.25)
time_elapsed4 = (time.perf_counter() - time_start)


np.mean(np.array([time_elapsed1,time_elapsed2,time_elapsed3,time_elapsed4]))
np.std(np.array([time_elapsed1,time_elapsed2,time_elapsed3,time_elapsed4]))



scipy.stats.pearsonr(analitics[:,1],X[:,1])
scipy.stats.pearsonr(analitics2[:,1],X2[:,1])
scipy.stats.pearsonr(analitics3[:,1],X3[:,1])
scipy.stats.pearsonr(analitics4[:,1],X4[:,1])

rms_data(analitics[:,1],X[:,1])
rms_data(analitics2[:,1],X2[:,1])
rms_data(analitics3[:,1],X3[:,1])
rms_data(analitics4[:,1],X4[:,1])

#####################################
fig, axs = plt.subplots(2, 2)
fig.suptitle('Solución usando Runge-Kutta (Cuarto orden)')
axs[0, 0].plot(X[:,0],X[:,1])
axs[0, 0].set_title('r=0.1')
axs[0, 1].plot(X2[:,0],X2[:,1], 'tab:orange')
axs[0, 1].set_title('r=0.5')
axs[1, 0].plot(X3[:,0],X3[:,1], 'tab:green')
axs[1, 0].set_title('r=1')
axs[1, 1].plot(X4[:,0],X4[:,1], 'tab:red')
axs[1, 1].set_title('r=5')


for ax in axs.flat:
    ax.set(xlabel='Tiempo (Minutos', ylabel='Individuos')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    
#####################################
    
fig, axs = plt.subplots(2, 2)
fig.suptitle('Solución analitica')
axs[0, 0].plot(analitics[:,0],analitics[:,1])
axs[0, 0].set_title('r=0.1')
axs[0, 1].plot(analitics2[:,0],analitics2[:,1], 'tab:orange')
axs[0, 1].set_title('r=0.5')
axs[1, 0].plot(analitics3[:,0],analitics3[:,1], 'tab:green')
axs[1, 0].set_title('r=1')
axs[1, 1].plot(analitics4[:,0],analitics4[:,1], 'tab:red')
axs[1, 1].set_title('r=5')


for ax in axs.flat:
    ax.set(xlabel='Tiempo (Minutos)', ylabel='Individuos')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    
    
#####################################
    
fig, axs = plt.subplots(4, 2,figsize=(8,10))
fig.suptitle('Solución usando Runge-Kutta (Segundo orden)')
axs[0, 0].plot(X[:,0],X[:,1])
axs[0, 0].set_title('r=0.1 Sol analitica')
axs[0, 1].plot(analitics[:,0],analitics[:,1])
axs[0, 1].set_title('r=0.1 Runge-Kutta (Segundo orden)')

axs[1, 0].plot(X2[:,0],X2[:,1], 'tab:orange')
axs[1, 0].set_title('r=0.5 Sol analitica')
axs[1, 1].plot(analitics2[:,0],analitics2[:,1], 'tab:orange')
axs[1, 1].set_title('r=0.5 Runge-Kutta (Segundo orden)')

axs[2, 0].plot(X3[:,0],X3[:,1], 'tab:green')
axs[2, 0].set_title('r=1 Sol analitica')
axs[2, 1].plot(analitics3[:,0],analitics3[:,1], 'tab:green')
axs[2, 1].set_title('r=1 Runge-Kutta (Segundo orden)')

axs[3, 0].plot(X4[:,0],X4[:,1], 'tab:red')
axs[3, 0].set_title('r=5 Sol analitica')
axs[3, 1].plot(analitics4[:,0],analitics4[:,1], 'tab:red')
axs[3, 1].set_title('r=5 Runge-Kutta (Segundo orden)')

for ax in axs.flat:
    ax.set(xlabel='Tiempo (Minutos)', ylabel='Individuos')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()    