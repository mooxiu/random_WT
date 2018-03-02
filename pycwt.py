import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from pandas import Series, DataFrame
import seaborn 
import TFwindow as TF
import wavelets 

# normalize the earthquake data and save to Series and treat them as function
f=np.loadtxt(r'C:\Users\monet\Dropbox\Research\20171214 test on python\ft\filtered.csv')
time=f[:len(f)-1,0]
real=f[:len(f)-1,1]
signal=Series(real, index= time)

# wavelet transform
# scale and time shift
sca=1
tshi=10.003
# time-frequency window
TC, TR=TF.Twindow(wavelets.WShannon, tshi, sca)


# Wavelet transform
dt= 0.01
i=np.arange(0, 119, dt) 
W1=np.sum(signal[i]* wavelets.WShannon(i, tshi, sca)*dt)
W2=np.sum(signal[i]* wavelets.WShannon(i, tshi+2*TR, sca)*dt)
