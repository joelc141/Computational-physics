# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 16:22:37 2018

@author: joelc
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 12:37:07 2018

@author: joelc
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import inv 
def f(g):
	def g(x):return g
	bohr=1
	nx=100
	nt=100
	t0=0
	tf=15
	tf=tf*2
	a=0.01*bohr
	b=15
	
	x=np.linspace(a,b,nx)
	t=np.linspace(t0,tf,nt)
	
	
	bc=0
	dx=(b-a)/nx
	dt=(tf-t0)/nt*1j
	numprotons=1
	#-(1.6)*10**(-19)*numprotons*((1.6)*10**(-19))*(8.99*10**28)
	# -(1.6)*10**(-19)*numprotons*((1.6)*10**(-19))*(8.99*10**9)*(10**18)
	k=-1
	#1.0*10**(-34)*10**18
	hb=1
	m=1
	#9.1*10**-31
	
	
	
	C=-dt/dx**2*(-hb/(2*m))
	alpha=C
	stability=alpha*dt/dx-1/2 #must be less than zero! or won't work
	
	
	#imaginary numbers
	
	def f(xx): return k/xx #bounding potential
	 #initial condition
	
	
	U=np.zeros( (nx,nt) )*1j
	PSI=np.zeros( (nx,nt) )*1j
	U[:,0]=g(x)
	
		
	v1=np.zeros(nx)*1j
	v2=np.zeros(nx-1)*1j
	v3=np.zeros(nx-1)*1j
	
	t0=hb**2/(2*m*dx**2)
	
	
	for i in range (0,nx):
		v1[i]=1/hb*(dt*(2*t0+f(x[i])))+1
	for i in range (0,nx-1):
		v2[i]=-1/hb*t0*dt
		v3[i]=-1/hb*t0*dt
		
	AA=np.zeros( (nx,nx) )*1j
	A1=np.diag(v1)
	A2=np.diag(v2,1)
	A3=np.diag(v3,-1)
	AA=A1+A2+A3
	
	
	for i in range (0,nt-1):
		U[:,i+1]=np.linalg.solve(AA,U[:,i])
	
	PSI=U
	
	PSI1=np.abs(PSI.conj()*PSI) #normalizing
	PSI=PSI1
	for j in range(0,nt):
		A=0
		for i in range(0,nx):
			A=A+PSI1[i,j]
		A=1/(dx*A)
		for k in range(0,nx):
			PSI1[k,j]=PSI1[k,j]*A
			return PSI1
