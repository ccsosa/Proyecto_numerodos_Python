# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 08:49:25 2019

@author: cami_
"""

# EDO. Método de Taylor 3 términos 
# estima la solucion para muestras espaciadas h en eje x
# valores iniciales x0,y0
# entrega arreglo [[x,y]]
########################################
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import sys
import os
import time

sys.path.append(os.path.abspath("E:/JAVERIANA/COMPUTACION/code_proyecto2/RMS.py"))
from RMS import*
#from file import rmsd, rms_data
#from RMS import*

########################################

def edo_taylor3t(d1y,d2y,x0,y0,h,xmax,r):
    intervals = int(xmax/h)
    print("intervals:"+str(intervals))
    estimado = np.zeros(shape=(intervals,2),dtype=float)
    #tamano = muestras + 1
    #estimado = np.zeros(shape=(tamano,2),dtype=float)
    # incluye el punto [x0,y0]
    estimado[0] = [x0,y0]
    x = x0
    y = y0
    for i in range(1,intervals,1):
        y = y + h*d1y(x,y,r) + ((h**2)/2)*d2y(x,y,r)
        x = x+h
        estimado[i] = [x,y]
        
        print(y)
    return(estimado)
    
    
 #  2946519
    
    # PROGRAMA PRUEBA
# Ref Rodriguez 9.1.1 p335 ejemplo.
# prueba y'-y-x+(x**2)-1 =0, y(0)=1
########################################
# INGRESO.
# d1y = y' = f, d2y = y'' = f'

#d1y = lambda x0,y0:(y0*x0)
#d2y = lambda x0,y0: (y0*y0)*x0


d1y = lambda x0,y0,r,:r*y0
d2y = lambda x0,y0,r: (r*r)*y0
    
#d1y = lambda x0,y0,: r*(y0*x0)
#d2y = lambda x0,y0,r: (r*r)*(y0*x0)

r=0.1
x0 = 0.00
y0 = 1E6
h = 0.25
xmax=200
########################################
# PROCEDIMIENTO
puntos = edo_taylor3t(d1y,d2y,x0,y0,h,xmax,r)
xi = puntos[:,0]
yi = puntos[:,1]

# SALIDA
print('estimado[xi,yi]')
print(puntos)


########################################
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
    
r=0.1   
analitics = growth(t0=0,N0=1E6,xmax=200,h=0.25,r=r)

# ERROR vs solución conocida
errores = analitics[:,1] - yi
errormax = np.max(np.abs(errores))

# SALIDA
print('Error máximo estimado: ',errormax)
print('entre puntos: ')
print(errores)

"""
# GRAFICA [a,b+2*h]
a = x0
b = h*muestras+2*h
muestreo = 10*muestras+2
xis = np.linspace(a,b,muestreo)
yis = y_sol(xis)
"""
xis = analitics[:,0]
yis = analitics[:,1]
# Gráfica



r=0.1   
X = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
analitics = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)
r=0.5
X2 = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
analitics2 = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)
r=1
X3 = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
analitics3 = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)
r=5
X4 = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
analitics4 = growth(t0=0,N0=1E6,xmax=20,h=0.25,r=r)


#Learn more or give us feedback
"""
N_0 = input('Give initial population size N_0: ')
r   = input('Give net growth rate r: ')
dt  = input('Give time step size: ')
N_t = input('Give number of steps: ')
"""
time_start = time.perf_counter()
r=0.1   
ff_list1 = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
time_elapsed1 = (time.perf_counter() - time_start)

time_start = time.perf_counter()
r=0.5   
ff_list2 = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
time_elapsed2 = (time.perf_counter() - time_start)

time_start = time.perf_counter()
r=1  
ff_list3 = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
time_elapsed3 = (time.perf_counter() - time_start)

r=5   
ff_list4 = edo_taylor3t(d1y,d2y,x0=0,y0=1E6,xmax=20,h=0.25,r=r)
time_elapsed4 = (time.perf_counter() - time_start)


np.mean(np.array([time_elapsed1,time_elapsed2,time_elapsed3,time_elapsed4]))
np.std(np.array([time_elapsed1,time_elapsed2,time_elapsed3,time_elapsed4]))




scipy.stats.pearsonr(analitics[:,1],X[:,1])
scipy.stats.pearsonr(analitics2[:,1],X2[:,1])
scipy.stats.pearsonr(analitics3[:,1],X3[:,1])
scipy.stats.pearsonr(analitics4[:,1],X4[:,1])

rms_data(np.array(analitics[:,1]),np.array(X[:,1]))
rms_data(np.array(analitics2[:,1]),np.array(X2[:,1]))
rms_data(np.array(analitics3[:,1]),np.array(X3[:,1]))
rms_data(np.array(analitics4[:,1]),np.array(X4[:,1]))


plt.scatter(analitics[:,1].T,X[:,1].T,c='red')
plt.plot(analitics[:,1].T,X[:,1].T,
         markerfacecolor='blue', markersize=12, color='skyblue', linewidth=2)
plt.plot(analitics2[:,1].T,X2[:,1].T,color="orange")
plt.plot(analitics3[:,1].T,X3[:,1].T,color="green")
plt.plot(analitics4[:,1].T,X4[:,1].T,color="red")


#####################################
fig, axs = plt.subplots(2, 2)
fig.suptitle('Solución usando Series de Taylor')
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
fig.suptitle('Solución usando Series de Taylor')
axs[0, 0].plot(X[:,0],X[:,1])
axs[0, 0].set_title('r=0.1 Sol analitica')
axs[0, 1].plot(analitics[:,0],analitics[:,1])
axs[0, 1].set_title('r=0.1 Serie de Taylor')

axs[1, 0].plot(X2[:,0],X2[:,1], 'tab:orange')
axs[1, 0].set_title('r=0.5 Sol analitica')
axs[1, 1].plot(analitics2[:,0],analitics2[:,1], 'tab:orange')
axs[1, 1].set_title('r=0.5 Serie de Taylor')

axs[2, 0].plot(X3[:,0],X3[:,1], 'tab:green')
axs[2, 0].set_title('r=1 Sol analitica')
axs[2, 1].plot(analitics3[:,0],analitics3[:,1], 'tab:green')
axs[2, 1].set_title('r=1 Serie de Taylor')

axs[3, 0].plot(X4[:,0],X4[:,1], 'tab:red')
axs[3, 0].set_title('r=5 Sol analitica')
axs[3, 1].plot(analitics4[:,0],analitics4[:,1], 'tab:red')
axs[3, 1].set_title('r=5 Serie de Taylor')

for ax in axs.flat:
    ax.set(xlabel='Tiempo (Minutos)', ylabel='Individuos')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
    
    

   