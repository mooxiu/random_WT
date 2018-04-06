# random_WT
This will help to change a certain time and frequency of a signal using wavelet transform

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




[1]:Fundmentals of Wavelets, p118
