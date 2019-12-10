# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 11:21:38 2019

@author: cami_
"""
import numpy as np
import matplotlib.pyplot as plt
import time

time_start = time.perf_counter()                     
##% Datos de entrada. Condiciones iniciales

Lx= 1                #% Ancho de la placa (m)
Ly= 3                #% Largo de la placa (m)
Nx=40                #% Número de puntos en el eje x
Ny=80                #% Número de puntos en el eje 
T_inf   = 30         #% Temperatra en el lado inferior ("Dirichlet Conditions" )
T_sup   = 80#0,20,50,80        #% Temperatra en en lado superior ("Dirichlet Conditions" )
T_der  = 30          #% Temperatra en en lado derecho ("Dirichlet Conditions" )
T_izq  = 30          #% Temperatra en en lado izquierdo ("Dirichlet Conditions" ) 


##% 

k=0                                              #% Contador de iteraciones
                                        #% delta y=h

#BOUNDARIES
Temp_v=np.zeros((Nx+2,Ny+2),dtype="float")
Temp_v2=np.zeros((Nx+2,Ny+2),dtype="float")
Temp_v[:,0]=T_izq;             Temp_v2[:,0]=T_izq
Temp_v[:,(Ny+1)-1]=T_der;      Temp_v2[:,(Ny+1)-1]=T_der
Temp_v[:,(Ny+2)-1]=T_der;          Temp_v2[:,(Ny+2)-1]=T_der            
Temp_v[(Nx+1)-1,:]=T_inf;          Temp_v2[(Nx+1)-1,:]=T_inf
Temp_v[(Nx+2)-1,:]=T_inf;          Temp_v2[(Nx+2)-1,:]=T_inf              
Temp_v[0,:]=T_sup ;            Temp_v2[0,:]=T_sup

##% 3- Steady-State section
contfig = 0;k=0
for x in range(0,100,1):
    print("k=",k)  
    #print("err_SS_max=",err_SS_max[k])    
    #print("Tolerancia_ss=",tolerencia_ss) 
    #for i=2:Nx                                                    
    for i in range(1,Nx-1,1):
        #for j=2:Ny
        for j in range(1,Ny-1,1):
            Temp_v2[i,j]=0.25*(Temp_v[i+1,j]+Temp_v[i,j+1]+Temp_v[i-1,j]+Temp_v[i,j-1])
        
    
    k=k+1                                                        
    contfig = contfig + 1
    #    err_SS_max(k)=abs(max(max(Temp_v2-Temp_v)))                         
    #err_SS_min(k)=abs(min(min(Temp_v2-Temp_v)))                        
    x=Temp_v2-Temp_v
    x_max = x;x_max=x_max.max(axis=0)
    print("x_max=",x_max)
    x_min = x;x_min=x_min.min(axis=0)

    Temp_v=Temp_v2   


    if contfig == 100:
        #fig = plt.figure()
        plt.imshow(Temp_v,cmap="jet")
        plt.colorbar()
        plt.title('Temp 80 Celsius')
        plt.grid(True)
        plt.show()
        #plt.imshow(Temp_v, interpolation ="nearest",norm=T)

 

time_elapsed = (time.perf_counter() - time_start)
time_elapsed