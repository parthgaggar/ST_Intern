from numpy import *
from scipy import signal
from pylab import *
from numpy import fft

N=500
fs = 500            
fc = 30                   
hd = zeros(N+1)
samples=arange(0,1.0,1.0/fs)
f_max = 50
step_size = 0.1
x = [1 + sin(2*pi*f*samples) for f in arange(0,f_max,step_size)]
phase_difference = zeros(f_max/step_size)
amplitude_gain = zeros(f_max/step_size)
wc = 2*pi*fc/fs
for n in arange(N+1):
    if n!=N/2:
        hd[n] = 1.0*sin(wc*(n-N/2))/(pi*(n-N/2))
    else:
        hd[n] = 1.0*wc/pi
    hd[n] = hd[n]*(0.54-0.46*cos(2.0*pi*n/N))
    
#plot (hd)
#show()
hd_fft = fft.rfft(hd)

plot (20*log10(abs(hd_fft)))
show()
#plot (angle(hd_fft,deg=True))
#show()
count = 0
for each_x in x:
    phase_diff = 0
    y= convolve(each_x,hd)
    y_fft = fft.fft(y)
    output = [y[b] for b in arange((fs+N-1)/2-fs/2,(fs+N-1)/2+fs/2,1)]
    output_fft = fft.fft(output)
    x_fft = fft.fft(each_x)
    #phase_diff = angle(output_fft,deg=True)-angle(x_fft,deg=True)
    #phase_difference[count] = mean(phase_diff)
    amplitude_gain[count] = 1.0*max(output)/max(each_x)
    count+=1
#plot (20*log10(amplitude_gain))
#show()
#plot (abs(phase_difference))
#show()
