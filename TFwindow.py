'''
This Program is used for calculate the RMS Radius of a time-frequency Window.
Twindow function will return values of Time center and Time radius.
Fwindow function will return values of Frequency center and Frequency radius.
'''


import pywt
import numpy as np
import matplotlib.pyplot as plt

'''
def MeyerScaling(t):
	for i in t==0:
		if i:
			v=2/3+4/(3*np.pi)
		else:
			v=(np.sin(2*np.pi*t/3)+4*t*np.cos(4*np.pi*t/3)/3)/(np.pi*t-16*np.pi*t**3/9)
	return v

def meyeraux(x):
	return 35*x**4-85*t**5+70*x**6-20*x**7
'''

# definition of functions
'''
def psi1(t):
	above=4/(3*np.pi)*(t-0.5)*np.cos(2*np.pi/3*(t-0.5))-1/np.pi*np.sin(4*np.pi/3*(t-0.5))
	bellow=(t-0.5)-16/9*(t-0.5)**3
	return above/bellow

def psi2(t):
	above=8/(3*np.pi)*(t-0.5)*np.cos(8*np.pi/3*(t-0.5))+1/np.pi*np.sin(4*np.pi/3*(t-0.5))
	bellow=(t-0.5)-64/9*(t-0.5)**3
	return above/bellow

def Meyer(t):
	psi=psi1(t)+psi2(t)
	return psi

def WindowMeyer(t,s):
	return Meyer(t/s)
'''
def Shannon(t):
	return np.sin(np.pi*t)/(np.pi*t)
# window function
# tshi=time shift, sca=scale
def WShannon(t, tshi, sca):
	return sca**(-0.5)*Shannon((t-tshi)/sca)




#Time radius and time center 
#dt is the time 

#time center and time radius
def Twindow(func, tshi, sca):
	dt=0.01
	n=30
	i=np.arange(-n,n,dt)    #time list
	norm2= np.sum(func(i, tshi, sca)**2*dt)
	norm= np.sqrt(norm2)
	
	temp1=np.sum(func(i, tshi, sca)**2*dt*i)
	TC= temp1/norm2
	
	temp2=np.sum((i-TC)**2*func(i, tshi, sca)**2*dt)
	TR= np.sqrt(temp2)/norm
	return TC, TR

# print(TR(WShannon, 9, 2), TR(WShannon, 0, 2))

# Frequency center and Frequency radius	
#Fourier Transform
def FT(func, tshi, sca):
	N=512
	dt=0.1
	
	t=np.arange(-N*dt,N*dt,dt)
	freq=np.linspace(0, 1/dt, 2*N)
	
	Amp=np.abs(np.fft.fft(func(t, tshi, sca)))
	
	freq=freq[:int(len(freq)/2)]
	Amp=Amp[:int(len(Amp)/2)]
	
	return freq, Amp

# 但是不同的dt的话，会有不同的频率的单位

def Fwindow(func, tshi, sca):
	freq, Amp= FT(func, tshi, sca)
	dw= freq[1]-freq[0]
	i= np.arange(0, len(Amp))
	
	fnorm2=np.sum(Amp[i]**2*dw)
	fnorm=np.sqrt(fnorm2)
	
	temp3= np.sum(freq[i]*Amp[i]**2*dw)
	FC= temp3/fnorm2

	temp4= np.sum((freq[i]-FC)**2*Amp[i]**2*dw)
	FR= np.sqrt(temp4)/fnorm

	return FC, FR
	


	

