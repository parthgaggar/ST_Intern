from numpy import *
from scipy import signal
from pylab import *
from numpy import fft

N=4096            #filter order
fs = 200         #sampling frequency=20 khz
fc = 25        #cutoff freq = 2.5 khz
w = ones(N+1) #rectangular window
hd = zeros(N+1)
wc = 2*pi*fc/fs
for n in arange(N+1):
    if n!=N/2:
        hd[n] = 1.0*sin(wc*(n-N/2))/(pi*(n-N/2))
    else:
        hd[n] = 1.0*wc/pi
h = [w[p]*hd[p] for p in arange(N+1)]
h_fft = fft.rfft(h)
plot (20*log10(abs(h_fft)))
show()
