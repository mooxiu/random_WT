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

This file is used to calculate the time and frequency window. For a given time shift and scale, we have to decide the center time $t^{*}$ and center scale $\omega ^{*}$,

[1]:Fundmentals of Wavelets, p118
