# In this file I will define some wavelet functions

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

def WMeyer(t, tshi, sca):
	return sca**(-0.5)*Meyer((t-tshi)/sca)


#define Shannon wavelet
def sinc(t):
	return np.sin(np.pi*t)/(np.pi*t)
def Shannon(t):
	return sinc(t/2)*np.cos(3*np.pi*t/2)
# window function
# tshi=time shift, sca=scale
def WShannon(t, tshi, sca):
	return sca**(-0.5)*Shannon((t-tshi)/sca)

