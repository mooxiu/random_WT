# random_WT
This will help to change a certain time and frequency of a signal using wavelet transform.



## Original Signal





## wavelets.py

This file is used to save wavelet functions. In my respository, the real part of Shannon Wavelet is used as an example. It can be used because real part of Shannon Wavelets are orthogonal[1]. According to [Wikipedia of Shannon Wavelet](https://en.wikipedia.org/wiki/Shannon_wavelet), the expression of the real part of Shannon Wavelets can be written as:

$\psi^{(Sha)}(t)= sinc(\dfrac{t}{2}) \cdot cos(\dfrac{3\pi t}{2})$

In which $sinc(t):= \dfrac{sin \pi t}{\pi t}$.

The window function for wavelet can be written as:

$\psi_{s,\tau}(t)= \dfrac{1}{\sqrt{a}} \psi(\dfrac{t-\tau}{s})$

In which s is scale, $\tau$ is time shift. In the program file, they are written as sca and tshi. And the window wavelet function is represented as WShannon. 

WShannon_sp is ???




## TFwindow.py

This file is used to calculate the time and frequency window. For a given time shift and scale, we have to decide the center time $t^{\*}$ and center scale $\omega ^{\*}$, the RMS radius of time $\Delta_{\phi}$ and RMS radius of scale $\Delta_{\hat{\phi}}$. According to [This slide](http://www.spcom.ecei.tohoku.ac.jp/~aito/wavelet/slide2.pdf), those parameters can be calculated by:

Time center: $t^{\*}:= \dfrac{1}{||\phi||^2} \int_{-\infty}^{\infty} t |\phi(t)|^2 dt$

RMS radius of time: $\Delta_{\phi}:= \dfrac{1}{||\phi||}\[\int_{-\infty}^{\infty}(t-t^{\*})^2|\phi(t)|^2 dt\]^{1/2}$

The center and radius for frequency is similar to time, they can be represented as:

Frequecny center: $\omega^{\*}:=\dfrac{1}{||\hat{\phi}||^2} \int_{-\infty}^{\infty} \omega |\hat{\phi}(\omega)|^2 d\omega$

RMS radius of frequency: $\Delta_{\hat{\phi}}:= \dfrac{1}{||\hat{\phi}||}\[\int_{-\infty}^{\infty} (\omega-\omega^{\*})^2|\hat{\phi}(\omega)|^2 d\omega\]$

In which, $||\phi||^2$ means the 2 norm of $\phi$, that is $||\phi||^2=\int_{-\infty}^{\infty}|\phi(t)|^2 dt$. 

This file are using the formulas above to calculate the important parameters of time-frequency windows.



## pycwt.py

This file can be divided into three parts, they are dealing with 1.wavelet transform, 2.phase calculation and randomization, 3. reconstruct of the signal.

### wavelet transform

Continuous wavelet transform can cover the whole time and frequency domian, which is unpractical. As we are only considering change a certain part of time and frequency, we can only integrate in a certain range to find out the wavelet coefficient of certain time and frequency.

The wavelet coefficient for all time range and frequency range:

$W_{f}(b,a)=\dfrac{1}{\sqrt{a}}\int_{-\infty}^{\infty} f(t) \bar{\psi}\dfrac{t-b}{a}dt$

If we have assign a certain time and frequency and want to figure out its wavelet coefficient and the wavelet coefficient which is consecutive to it. We have to notice, that the two wavelet have the same frequency but different time shift. However, same frequency means that they have the same RMS time radius. Thus we assume that the distance of center time of two consecutive wavelets is twice the RMS time radius.

We assume to change the signal in scale of 1 and time shift of 10.003( this is to hope that in the later calculation we won't divide by 0).

$W_1$ and $W_2$ are two consecutive wavelet coefficients:

```python
dt= 0.01
i=np.arange(0, 119, dt) 
W1=np.sum(signal[i]* wavelets.WShannon(i, tshi, sca)*dt)
W2=np.sum(signal[i]* wavelets.WShannon(i, tshi+2*TR, sca)*dt)
```

To find out the corresponding signal of these two parts, we need to do inverse wavelet transform of them. The inverse wavelet transform is :

$f(t)=\dfrac{1}{C_{\psi}}\int_{-\infty}^{\infty} db \int_{-\infty}^{\infty} \dfrac{1}{a^2}[W_{\psi}f(b,a)] \psi_{b,a}(t) da$

We need to calculate the $C_{\psi}$ which is  a constant first. The presentation of $C_{\psi}$ is like:

$C_{\psi}= \int_{-\infty}^{\infty} \dfrac{|\hat{\psi}(\omega)|^2}{|\omega|} d\omega$

```\python
freq, Amp= TF.FT(wavelets.WShannon, tshi, sca)
dw=freq[1]-freq[0]
i=np.arange(0, len(Amp))
Cpsi=np.sum(Amp[i]**2/freq[i]*dw)
```

This is the Fourier Transform of the wavelet:

![img](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAd8AAAFJCAYAAADaPycGAAAABHNCSVQICAgIfAhkiAAAAAlwSFlz%0AAAALEgAACxIB0t1+/AAAHoVJREFUeJzt3X+QldWd5/HPc+/t393QDVxM1OCCWbKaycQxtRprRStl%0AEJNZE60lJZBtkyL7x7pMGVKJIpTpmJAiUP5hMmwBxko2O+hGnWhlSO1uMquhYiQO6ziio4MkSxwi%0AyJAGGvsn3ffe59k/bt+nG20auPc5zznn9vtVZXXT3h9fDg2f/p7zPOcEURRFAgAAqcnYLgAAgJmG%0A8AUAIGWELwAAKSN8AQBIGeELAEDKCF8AAFKWS+NNensHEn/Nrq5W9fUNJ/66MwljWDvGsHaMYTIY%0Ax9olPYb5fMdZ/5+3nW8ul7VdgvcYw9oxhrVjDJPBONYuzTH0NnwBAPAV4QsAQMoIXwAAUkb4AgCQ%0AMsIXAICUEb4AAKSM8AUAIGWELwAAKSN8AQBIGeELAEDKCF9D3jzar9++dcp2GQAABxG+hmz873+v%0AzY/9g+0yAAAOInwBAEgZ4WtYsRTaLgEA4BjC17DTYyXbJQAAHEP4GnZ6tGi7hPeIokhRFNkuAwBm%0ALMLXMBc73w2P7NVDf/2K7TIAYMbK2S6g3rkWvlEU6djJYR07OWy7FACYseh8DTs95ta087CD0+AA%0AMNMQvoYNjhSc2mxjcLhguwQAmPGYdjbsv/3vN1QohvrPn/2wrrniItvlaIDwBQDrzqvzfeWVV9Td%0A3S1JOnTokFauXKlVq1bpG9/4hsKQ+1inUyiWx2fH37yuH/6v/ZarkQaGx2yXAAAz3jk730ceeUS7%0Adu1SS0uLJOk73/mO1q5dq2uvvVY9PT169tlntXTpUuOF1oPnXz2qf/cn79PfvviWrvxXc9TR2qBr%0ArrhIURSpFEYKAimbSX4l4PlXj2r4dEH/dKhP+dktib8+AODCnDN8FyxYoK1bt+ree++VJL3++uu6%0A5pprJEk33HCD9uzZQ/hegC3/42VJ0su/Oy5JOnikX3/3T/+iweGCFl0ySzf86cU6+Ha/CsVQs9oa%0AtGB+h04XSjo9WlRzY1aNDVk15DLKBIGymUDDo0UViqHCKFIuWw7uNw71aWS0qJbmnN7X1aqfPv+m%0Atd8vAOC9zhm+y5Yt0+HDh+NfR1GkIAgkSW1tbRoYGDBX3Qzwf/7+rfjzg0f6dfBIv8VqAABpuOAL%0ArjKTpkWHhoY0a9ascz6nq6tVuVz2Qt/qnPL5jsRfM20fXjRXr//+hCTp0vntumpxXqdHS3rj0El9%0A4KIOXX7pbBWKoU6+c1ptLQ3Kd7UoGwQqhZHGiqECSf1DYyqGoeZ0NGtgeEz7//mk3jo2oMULujR8%0Auqj9/3zyrO9fD2NoG2NYO8YwGYxj7dIawwsO3yuvvFJ79+7Vtddeq+eee04f//jHz/mcvr7kN3TI%0A5zvU2+tu192Yy2iseObFaJ+6doFuuOpi/c8XDun/7j+mh/7ierU05dQ3MKpTg6Na+P5z/yBzPv78%0A2gUKo0iZIFAUReobGNV9D7+gGz56sV5/86SO9Y1IksIw0okTg4m850zl+vehDxjDZDCOtUt6DKcL%0A8gu+umfdunXaunWr7rjjDhUKBS1btqym4upVc1P555rrP/L++Gt/tjivi7pa9YVbPqS/vHuJWsYf%0A09XRlFjwVmTGlwaCINCcWc36r2tv0KpPLtZf/Ic/jR9TCtnfGQBsOK/O99JLL9WTTz4pSVq4cKEe%0AffRRo0XVhSjS++a0avWfX6E3j/bryPEhzZvdLKl8RXM25e1NGhvK0/6XzGvTRxbN1T/+/oRKHHcI%0AAFawyYYhkaTx5lP3/cerdWpwTJ3tTVZrqshmyoUV6XwBwArC15DJJ/a1NTeorbnBXjHvks2Ww5fO%0AFwDsYG9ngyq3ZLmm0vmy5gsAdhC+hpTvh7ZdxdTiaWc6XwCwgvA1JIokR7M33sKyVKLzBQAbCF9D%0AyrHmZvxW1nzpfAHADsLXGHennTPj084ha74AYAXha4jb0850vgBgE+FrSCQ5m765ypovnS8AWEH4%0AmhJJgaPpy5ovANhF+BoSyd1558q+z3S+AGAH4WuKu9nLDlcAYBnha8jkvZ1dM3HBFZ0vANhA+BoS%0AOXzF1cQmG3S+AGAD4WuMu/f5srczANhF+Bri9H2+8Zov4QsANhC+JjmavhPn+TLtDAA2EL6GRC7f%0A58vBCgBgFeFriMv3+U6s+dL5AoANhK8hUeTu4E7scEXnCwA2uJoP9cHRy53jzpdbjQDACsLXgChy%0A+TTfSWu+3GoEAFYQvgZUIs3Rxjc+z5eDFQDADsLXBMcbyvg+XzpfALCC8DUgGk/fwNHWN0fnCwBW%0AEb4GRI43lJU135DOFwCsIHwNcrTxnbTmS/gCgA2ErwGVztfR7OVWIwCwjPA1opK+bsbvxCYbhC8A%0A2ED4GuBN58uaLwBYQfgaEEeao+k7Me1M+AKADYSvCXHn62b6Vm6BCl2/LBsA6hTha8DEfb6WCzkL%0AR8sCgBmD8DXA+YZyPH0j5wsFgPpE+BqUcbT1rUyHE70AYAfha4DrHaWjPxMAwIxB+Brg+qlGMbd/%0ARgCAukX4GuB44xvzpU4AqDeEr0GunmpUqSui9QUAKwhfAyprvm5G7yRkLwBYQfga4PoOV5WGnOwF%0AADsIXxMc39u5UpfrV2UDQL0ifA2Y6Hwdjd+A+3wBwCbC1wTH13zjukhfALAiV82TCoWC7rvvPh05%0AckSZTEYbN27U5ZdfnnRt3nL+Pt94zZf0BQAbqup8f/WrX6lYLOrxxx/XmjVr9N3vfjfpurzm+lLq%0AxJqv1TIAYMaqKnwXLlyoUqmkMAw1ODioXK6qBrruuX6fLwDAjqpSs7W1VUeOHNGnPvUp9fX1aceO%0AHdM+vqurVblctqoCp5PPdyT+mkkIGsrD2tzU4GSNzcNjkspXO7tYn28Yw9oxhslgHGuX1hhWFb4/%0A+tGPdP311+urX/2qjh49qi984Qv62c9+pqampikf39c3XFORU8nnO9TbO5D46ybhZP9pSdLoWMHJ%0AGodOFySVp51drM8nLn8f+oIxTAbjWLukx3C6IK8qfGfNmqWGhgZJ0uzZs1UsFlUqlaqrrg5F8X2+%0Abk7vulkVAMwcVYXvF7/4RW3YsEGrVq1SoVDQV77yFbW2tiZdm7fi7SWdTTlnCwOAGaGq8G1ra9P3%0Avve9pGupG/GtRlarODeudgYAO9hkwwB/9nYmfQHABsLXhHiHK0fTdxydLwDYQfga4EvnCwCwg/A1%0AwflTjcYPVqD1BQArCF8D/NnbGQBgA+FrwERH6Wb6cqoRANhF+BrkaucbX+3MtDMAWEH4GhA5vuZb%0AqYzoBQA7CF8DJq52djN+47JIXwCwgvA1Ib7P121ssgEAdhC+Bvhyny9LvgBgB+FrguNrvq7vvAUA%0A9Y7wNWDiPl9HQ47OFwCsInwNiBxf85243or0BQAbCF8DHN9jAwBgGeFrkKtrq5XpcKadAcAOwteA%0AynSuq0u+AAC7CF8DfOgoA7G9JADYQvga5HTnG/jxQwIA1CPC14CJvZ3dTV+XawOAekf4GhC5vsuG%0Ayl05084AYAfha4L72SuJcxUAwBbC1wDX93aWxtejSV8AsILwNcGDNV8pYIcrALCE8DXAh/t8A652%0ABgBrCF8DfAg1Zp0BwB7C1yCXO9/xXTZsVwEAMxLha8DEqUbupm+ggM4XACwhfA2IN9lwN3vZ4QoA%0ALCJ8DSDTAADTIXxNqEw7O9z6BhI/JQCAJYSvAZVMczd6x281In0BwArC1wAfdriSAtZ8AcASwtcE%0AD/Z2drk2AKh3hK8BE6cauRtxnGoEAPYQviZ40PlKXG8FALYQvgbEF1w5nL5BwJovANhC+BrgT6h5%0AUygA1BXC1wgP7vNlhysAsIbwNSDyYM2XcxUAwB7C1wAv7vMNOFQQAGwhfE2IO19305fOFwDsIXwN%0AiOI1X8uFTCeg7wUAWwhfA3xZ8yV9AcCOXLVPfPjhh/XLX/5ShUJBK1eu1Oc+97kk6/Ja5MONvgAA%0Aa6oK37179+rll1/Wj3/8Y42MjOiHP/xh0nV5bnza2XIV0wmCgFONAMCSqsL3+eef1+LFi7VmzRoN%0ADg7q3nvvTbour0U+nCkoLrgCAFuqCt++vj69/fbb2rFjhw4fPqy77rpLP//5z8+6qURXV6tyuWxN%0AhU4ln+9I/DWTMOvtAUlSR3uzszVmsxlFcncMfcIY1o4xTAbjWLu0xrCq8O3s7NSiRYvU2NioRYsW%0AqampSSdPntTcuXOnfHxf33BNRU4ln+9Qb+9A4q+bhHf6RyRJQ0OjztYYhaEUZJytzxcufx/6gjFM%0ABuNYu6THcLogr+pq54997GP69a9/rSiKdOzYMY2MjKizs7PqAuuNH0f1Baz4AoAlVXW+n/jEJ/Ti%0Aiy9q+fLliqJIPT09ymaTn1b2ncsXO7O3MwDYU/WtRlxkdXY+3OcrifQFAEvYZMOA+BYeh1vfgB2u%0AAMAawtcEDzrfQAGNLwBYQvga4MepRhK9LwDYQfia4EXny5IvANhC+BowcaqRw/EbcKsRANhC+Brg%0AQ0fJrDMA2EP4GhC5f7Hz+NXOpC8A2ED4GhBPOzu96gsAsIXwNcGDzlfyY3ocAOoR4WuAD5kWBNzn%0ACwC2EL4Gudz5lksjfQHABsLXgMqpRk6v+XKwAgBYQ/ga4MMOV4HoewHAFsLXBA92uKL1BQB7CF8D%0AvOh8OdUIAKwhfE3wYM2XvZ0BwB7C14BKprl8tTP7SwKAPYSvAT50lJznCwD2EL4GuX2qEX0vANhC%0A+BowcZ+vu1jzBQB7CF8DfMi0gPQFAGsIXxPigxXc7n2JXgCwg/A1IJ52djh7A/bYAABrCF8D4luN%0ArFYBAHAV4WtA5EH6cqoRANhD+Brk8g5XbO0MAPYQvgZEHpysEHDBFQBYQ/ia4H72cqMvAFhE+Brg%0Aw97OnOcLAPYQvgZEHlxxReMLAPYQvga53Pm6XRwA1DfC14DIgzXfSm0R7S8ApI7wNSCOM4fTt9L4%0AEr0AkD7C14T4VCOH07eC9AWA1BG+BvjR+Vb2uCJ9ASBthK8JHqz5VrDkCwDpI3wNmLjP1934dbg0%0AAKh7hK8B8ZGCluuYzsTVzlbLAIAZifA1IL7VyOX0BQBYQ/ga4EUzGf9k4EW1AFBXCF8jxqedHW59%0AmXYGAHsIXwO8CDQ22QAAawhfgxxufCc2ACF9ASB1hK8BE3s7u5u+E9tLkr4AkLaawvfEiRO68cYb%0AdfDgwaTqqQuRR7tseDFFDgB1purwLRQK6unpUXNzc5L11AcPstfl2gCg3lUdvlu2bNGKFSs0f/78%0AJOupC17t7UznCwCpy1XzpKefflpz5szRkiVL9P3vf/+cj+/qalUul63mraaVz3ck/ppJaGlplCTN%0A6WpztsampvIf/dx57WpvabBcjd9c/TP2CWOYDMaxdmmNYVXh+9RTTykIAr3wwgvav3+/1q1bp+3b%0Atyufz0/5+L6+4ZqKnEo+36He3oHEXzcJw8OjkqRTp4bV21LVEBs3NlaUJB0/PqCRZsK3Wi5/H/qC%0AMUwG41i7pMdwuiCvKhkee+yx+PPu7m498MADZw3emcinqVyfagWAesGtRgY5fZ+vy8UBQJ2reU50%0A586dSdRRV7y4z3f8Y0TrCwCpo/M1ID5S0N3sZXtJALCI8DWAQAMATIfwNSE+z9fd1jeujJ8UACB1%0AhK8Ble0l3Y1exXPiZC8ApI/wNcCLHa4qn3DBFQCkjvA1wYe9nbngCgCsIXwNmOh8XY7fMhpfAEgf%0A4WtC5P6ar8sXgwFAvSN8Dag0ky7nG5tsAIA9hK8BXuSZwz8YAEC9I3yNqOxw5W7CTXS+VssAgBmJ%0A8DUg8uBq50p1Edc7A0DqCF8DvLjPN259rZYBADMS4WuCB50v2QsA9hC+BkQerPlyqhEA2EP4GuDD%0Amm981jBXXAFA6ghfA+I8czl9AQDWEL5GVHa4cjh9mXYGAGsIXwPiaWf3s5f0BQALCF8DSmE50bIZ%0Ad9OXU40AwB7C14CwcrCCw+Ebb7LBBVcAkDrC14DQo84XAJA+wteASvhmHE449nYGAHsIXwPGs9fp%0A8K20vmQvAKSP8DWgsuabcXnaufIJrS8ApI7wNSCednZ5dLnaGQCscTkevOXTmi/pCwDpI3wNCKNI%0AQeD2wQqBWPMFAFsIXwPCMHK665U0Me3Mmi8ApI7wNSCMIqcvtpI48wEAbCJ8DQhDt690liROFAQA%0AewhfA0o+TDsDAKwhfA2IokjuN76OFwgAdYzwNaAURk7v6yxNPtWIeWcASBvha0AYRY6faDSBNV8A%0ASB/ha4APtxo5Xh4A1DXC14Aw8mDaOT7P13IhADADEb4G+ND5ijVfALCG8DUgjNy/z5e9nQHAHsLX%0AgDB0f4crTjUCAHsIXwPK0862q5heQPoCgDWErwFe7O3Mmi8AWEP4GuDFBVfjuNoZANKXq+ZJhUJB%0AGzZs0JEjRzQ2Nqa77rpLN910U9K1ecunzhcAkL6qwnfXrl3q7OzUgw8+qFOnTum2224jfCfx4lSj%0A+D5fWl8ASFtV4XvLLbdo2bJlksr/eGez2USL8lkUReXO1/HW0u3qAKC+VRW+bW1tkqTBwUHdfffd%0AWrt27bSP7+pqVS6XfEDn8x2Jv2atSmG5k2xuyjlZX0VrW6MkafbsVqfr9AHjVzvGMBmMY+3SGsOq%0AwleSjh49qjVr1mjVqlW69dZbp31sX99wtW9zVvl8h3p7BxJ/3VoVS6EkqVQsOVlfxcjwmCTp1Klh%0Ap+t0navfhz5hDJPBONYu6TGcLsirCt/jx49r9erV6unp0XXXXVd1YfWo0vl6c6qR7QIAYAaq6laj%0AHTt2qL+/X9u2bVN3d7e6u7t1+vTppGvzUjgevs6v+cY3+hK/AJC2qjrf+++/X/fff3/StdSFcDzM%0A3D/VqIzoBYD0sclGwnzpfNldEgDsIXwTNp69zq/5cqoRANhD+Cas0vm6Pu1c2eKKvZ0BIH2Eb8Im%0App0tF3IOdL4AYA/hm7DKBVeur/kGrPkCgDWEb8Liztf11nccdxoBQPoI34TFna/j4Rvf50vvCwCp%0AI3wT5kvnG0cv2QsAqSN8E1biPl8AwDkQvgmrdJKuh29A+gKANYRvwnzZXrKC+3wBIH2Eb8ImTjWy%0AXAgAwFlERMJ82duZi50BwB7CN2G+bC9J9gKAPYRvwnzZ4Sre25n0BYDUEb4Jq4SvN6ca0fsCQOoI%0A34T5Mu0c32lE9gJA6gjfhIVh+aPr085uVwcA9Y3wTdjEmq/lQs4hiM/zBQCkjfBNmC97O1dEzDsD%0AQOoI34T5c6qR7QoAYOYifBPmzSYb4x9pfAEgfYRvwkq+TDvHa76kLwCkjfBNmC+bbMTVkb0AkDrC%0AN2HxkYKujywnCgKANa5HhHe8mXauIH0BIHWEb8J8u+AKAJA+wjdhhWJ5i6uGrNtDG3DBFQBY43ZC%0AeKh/eEyS1NHWaLmS88OtRgCQPsI3Yf1D5fCd5Xj4Mu0MAPYQvgmrhO/sVrfDl1ONAMAewjdh/UNj%0AamrMqqkxa7uUaQVizRcAbCF8E/bO0Jj7Xa8m7e1M9gJA6gjfBIVRpIHhgvPrvZORvQCQPsI3QYMj%0ABYVR5EX4On4bMgDUNcI3Qf2D4xdb+RC+lTVfrrgCgNQRvgk6cnxIkpTvbLFcyXlgb2cAsIbwTdDv%0ADp+SJP3rS2dbruTcONUIAOwhfBP0u8PvqCGX0WXv67BdyrnR+QKANYRvQnpPjejwHwd1+cWzlHN8%0AX2dpYs2XXTYAIH3up4Qn/vbFtxRJWvLRi22XckGIXgBIH+GbgAN/6NPufziiebOb9W//zXzb5ZyX%0Aro4mSdIfjg1YrgQAZh7CtwZhFOk3rx3Vd//6VUnSf/r3V3ox5SxJH7xktubMatJLB3pVLIW2ywGA%0AGSVnuwDfFEuh3vrjoF5786T+7vV/0dETw2pqyOq/3P4nWvyBTtvlnbdMJtCSqy7V3zx3UH/1iwO6%0Ac9mHvPnBAQB8V1X4hmGoBx54QAcOHFBjY6O+/e1v67LLLku6tlSEYaRCMVShFKpQDDUyWtTw6aIG%0ATxc0NFL+7+TAqI6/c1rHT43o7RNDKpbKK6W5bKDrPnyRbluyyI97e99l5c0f0ssHjun5V4/qwB/6%0A9Ik/u1QfWTRH75/XpgxbYAGAMVWF7zPPPKOxsTE98cQT2rdvnzZv3qzt27cnXduUwijSX/38gN4Z%0AHtPoaFFhVN6lKRr/OPnXYRSd8XkYRnHIVv4rhed/yVFDLqNL5rVr0cWz9MFLZuujH5yn1mZ/Jw/a%0AWhq0btXV+smvDur5V4/qyd3/T0/ulhpzGeU7WzRvdrM6WhvV1pJTW3OD2ppzamzIqiGXUUMuo8bc%0A5M8zymYzygTlrjqbmfx5oEwmUCYY/zj+OQDMVFUlx0svvaQlS5ZIkq666iq99tpriRY1nUIx1D/+%0A/oT6BkbjrwUq/yMfBFIQTHzMBFImCOKvZYJADbmM2lsa4tBoyGbUMClEmhqyam9pOCNwujqax4Oo%0AQUGdhUZLU07dN39Ity9ZpJcO/FG/O/yODvcOqvfUSLxjlwmVP7NMZvymp8qfncr7TgcK4v2n4z/T%0A8c816fMzvl55blB5zfL3wLtfY9IWIxPv8e7i4k/f+4DJj21oyKpYLE35mmc8dqrXnPT1ya95Zi3v%0Arnjqfbmn+r6cso6zvOZZ66vChf4VaWzMaWysWNN7nvH+tTzX8t/vWt4+iXGs9fdf8+jV+AK1PL2x%0AIasv3faR2n8P56mq8B0cHFR7e3v862w2q2KxqFxu6pfr6mpVLpfc+bY/6lmmKIrOCFpUJ58vbwiS%0Al7RwwZz461EUaeh0UQNDYxoYLv83OFzQWKGksUJJo4VQhWJJo4WSxgqhxoolhWGkUilSGJU/lsLy%0AzEIYRvHHyuelMIw/j8pvqEhSFJbPGK7MZETl/xXPYEjl2Q1F5ceF0eTnjs98hFH8XMUzIJOeO/Gb%0A1PhD3v2lyq/e87UpHzvpAVO/VjTF1yZ99ZzPB2BaJpCWXbdQV6d0x0pV4dve3q6hoYmuKAzDswav%0AJPX1DVfzNtPK5zvU28ttMrU4nzHMSepqyamrJSfNTacun6T9fRhN8QODzgj09/7AcObz3/uks/1w%0AUV2BF/6UefPadfz44PjTa6vA5g8stb93bS8wd267jp8YrP7dLf+wV+shL7WWn8tkdNkHuhL9+1xp%0AbqZ8v2pe8Oqrr9bu3bv16U9/Wvv27dPixYurLg7A+QvOMWWcwMRf6pqbcmpqTG5mbKZqb23UyFCD%0A7TJwnqoK36VLl2rPnj1asWKFoijSpk2bkq4LAIC6VVX4ZjIZfetb30q6FgAAZgR2VQAAIGWELwAA%0AKSN8AQBIGeELAEDKCF8AAFJG+AIAkDLCFwCAlBG+AACkjPAFACBlQVTrbtYAAOCC0PkCAJAywhcA%0AgJQRvgAApIzwBQAgZYQvAAApI3wBAEiZV+EbhqF6enp0xx13qLu7W4cOHbJdkrdeeeUVdXd32y7D%0AW4VCQffcc49WrVql5cuX69lnn7VdkndKpZLWr1+vFStWaOXKlfrtb39ruyRvnThxQjfeeKMOHjxo%0AuxRv3X777eru7lZ3d7fWr19v/P1yxt8hQc8884zGxsb0xBNPaN++fdq8ebO2b99uuyzvPPLII9q1%0Aa5daWlpsl+KtXbt2qbOzUw8++KBOnTql2267TTfddJPtsryye/duSdLjjz+uvXv36qGHHuLvcxUK%0AhYJ6enrU3NxsuxRvjY6OKooi7dy5M7X39Krzfemll7RkyRJJ0lVXXaXXXnvNckV+WrBggbZu3Wq7%0ADK/dcsst+vKXvyxJiqJI2WzWckX++eQnP6mNGzdKkt5++23NmjXLckV+2rJli1asWKH58+fbLsVb%0Ab7zxhkZGRrR69Wrdeeed2rdvn/H39Cp8BwcH1d7eHv86m82qWCxarMhPy5YtUy7n1aSHc9ra2tTe%0A3q7BwUHdfffdWrt2re2SvJTL5bRu3Tpt3LhRt956q+1yvPP0009rzpw5cVOC6jQ3N+tLX/qSfvCD%0AH+ib3/ymvva1rxnPFq/Ct729XUNDQ/GvwzAkRGDN0aNHdeedd+qzn/0swVGDLVu26Be/+IW+/vWv%0Aa3h42HY5Xnnqqaf0m9/8Rt3d3dq/f7/WrVun3t5e22V5Z+HChfrMZz6jIAi0cOFCdXZ2Gh9Hr8L3%0A6quv1nPPPSdJ2rdvnxYvXmy5IsxUx48f1+rVq3XPPfdo+fLltsvx0k9/+lM9/PDDkqSWlhYFQaBM%0Axqt/kqx77LHH9Oijj2rnzp264oortGXLFuXzedtleecnP/mJNm/eLEk6duyYBgcHjY+jV23j0qVL%0AtWfPHq1YsUJRFGnTpk22S8IMtWPHDvX392vbtm3atm2bpPKFbFz0cv5uvvlmrV+/Xp///OdVLBa1%0AYcMGxg9WLF++XOvXr9fKlSsVBIE2bdpkfFaVU40AAEgZczwAAKSM8AUAIGWELwAAKSN8AQBIGeEL%0AAEDKCF8AAFJG+AIAkDLCFwCAlP1/vB8M0Ay0E+cAAAAASUVORK5CYII=)

Then we can rewrite the presentation of inverse wavelet transform as:

$f(t)=\dfrac{W_{\psi}f(b,a)}{C_{\psi}} \int_{-\infty}^{\infty} db \int_{-\infty}^{\infty} \dfrac{\psi_{b,a}(t)}{a^2}da$



It's obvious that the left part of the sentence $\dfrac{W_{\psi}f(b,a)}{C_{\psi}}$is a constant. And we are going to integrate $\dfrac{\psi_{b,a}(t)}{a^2}$ over time shift $b$ and scale $a$.

In the pyfile $\dfrac{\psi_{b,a}(t)}{a^2}$ is marked as a new function WShannon_sp.

We calculate different time $t$ for the same time shift and frequency range. After that, we will get an array of $\int \int \dfrac{\psi_{b,a}(t)}{a^2} db da$.

We can find the pick up part of the signal.

The phase is assumed to be: $\theta= arctan (W_!/W_2)$, after randomizing the phase we get a new phase $\theta'$. The new wavelet coefficients will be $W_{1new}=\sqrt{W_1^2+W_2^2} \cdot sin\theta'$ and $W_{2new}=\sqrt{W_1^2+W_2^2} \cdot cos\theta'$.

And for the same reason, we will get the picked part of the signal after randomization.



[1]:Fundmentals of Wavelets, p118 
