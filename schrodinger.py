# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 12:37:07 2018

@author: joelc
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import inv 

nx=35
nt=100
t0=0
tf=0.0000000001
a=0.5
b=100

x=np.linspace(a,b,nx)
t=np.linspace(t0,tf,nt)


bc=0
dx=(b-a)/nx
dt=(tf-t0)/nt
numprotons=1
#-(1.6)*10**(-19)*numprotons*((1.6)*10**(-19))*(8.99*10**28)
k=-(1.6)*10**(-19)*numprotons*((1.6)*10**(-19))*(8.99*10**9)*(10**18)
hb=1.0*10**(-34)*10**18
m=9.1*10**(-31)



C=-dt/dx**2*(-hb/(2*m))
alpha=C
stability=alpha*dt/dx-1/2 #must be less than zero! or won't work




def f(x): return -k/x #bounding potential
def g(x): return np.sin(x) #initial condition


U=np.zeros( (nx,nt) )
PSI=np.zeros( (nx,nt) )
for i in range (1,nx):
	U[i,0]=g(x[i])
	PSI[i,0]=g(x[i])
#for j in range (0,nt-1):
#	for i in range(1,nx-1):
#		U[i,j+1]=C*(U[i+1,j]-2*U[i,j]+U[i-1,j])+U[i,j]*(dt/hb*f(x[i])+1)
#	
	
v1=np.zeros(nx)
v2=np.zeros(nx)
v3=np.zeros(nx)

t0=hb**2/(2*m*dx**2)


for i in range (0,nx):
	v1[i]=2*t0+f(x[i])
	v2[i]=-t0
	v3[i]=-t0
AA=np.zeros( (nx,nx) )
A1=np.diag(v1)
A2=np.diag(v2,1)
A3=np.diag(v3,-1)

for i in range (0,nx):
	for j in range(0,nx):
		AA[i,j]=A1[i,j]+A2[i,j]+A3[i,j]
AA=inv(AA)

for i in range (0,nt-1):
	PSI[:,i+1]=np.dot(AA,PSI[:,i])


		
PSI=abs(U*U) #normalizing

for j in range(0,nt):
	A=0
	for i in range(0,nx):
		A=A+PSI[i,j]
	A=1/(dx*A)
	for k in range(0,nx):
		PSI[k,j]=PSI[k,j]*A

			   
	
