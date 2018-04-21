import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from pandas import Series, DataFrame
import seaborn 
import TFwindow as TF
import wavelets 
import scipy.integrate as integrate

# normalize the earthquake data and save to Series and treat them as function
f=np.loadtxt(r'C:\Users\monet\Dropbox\Research\20171214 test on python\ft\filtered.csv')
time=f[:,0]
real=f[:,1]
signal=Series(real, index= time)

# wavelet transform
# scale and time shift
sca=1
tshi=10.003
# time-frequency window
TC, TR=TF.Twindow(wavelets.WShannon, tshi, sca)
FC, FR=TF.Fwindow(wavelets.WShannon, tshi, sca)

# Wavelet transform
dt= 0.01
i=np.arange(0, 119, dt) 
W1=np.sum(signal[i]* wavelets.WShannon(i, tshi, sca)*dt)
W2=np.sum(signal[i]* wavelets.WShannon(i, tshi+2*TR, sca)*dt)

# Calculation of Cpsi 
freq, Amp= TF.FT(wavelets.WShannon, tshi, sca)
dw=freq[1]-freq[0]
i=np.arange(0, len(Amp))
Cpsi=np.sum(Amp[i]**2/freq[i]*dw)

# reverse
inner=[]
for t in np.arange(0, 119, 0.1):
	def WShannon_sp(tshi, sca):
		return wavelets.WShannon(t, tshi, sca)/sca**2
	inner.append(integrate.nquad(WShannon_sp, [[0,119],[FC-FR,FC+FR]]))
inner_array=np.array(inner)

W1_wave=(W1/Cpsi)*inner_array
W2_wave=(W2/Cpsi)*inner_array
W_wave=W1_wave+W2_wave


# theta
theta= np.arctan(W1/W2)

# randomize theta
theta1= theta*np.random.rand()

# new W1 and new W2
W1new= np.sin(theta1)* np.sqrt(W1**2+W2**2)
W2new= np.cos(theta1)* np.sqrt(W1**2+W2**2)

# Reverse wavelet transform
W1new_wave=(W1new/Cpsi)*inner_array
W2new_wave=(W2new/Cpsi)*inner_array
Wnew_wave=W1new_wave+W2new_wave

#difference from the original signal
dWave=Wnew_wave-W_wave
# new signal
newsignal=Series(dWave+real, index=time)
plt.plot(newsignal)
plt.show()