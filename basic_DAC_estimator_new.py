from numpy import *
from pylab import *
from scipy import signal
from numpy import fft
from random import randint
import numpy as np

no_of_bits = 14  # Number of Bits
base=2     #radix
Vm = 8                          # Max Voltage
alpha = 0.001                   #higher order coefficient
order = 2                       #order of non linearity
################Sample and Hold Circuit###############
dig_out_adc = zeros(no_of_bits)
no_of_inputs = 10000
#################ADC successive approximation#########
f=1000.0
sample_length = 1.0
sample_duration = 0.0001
t = np.arange(0 ,sample_length, sample_duration)
x = 1024.0 + 512.0 * (sin(2.0*pi*f*t/sample_length))
#x = [1024.0 for t in arange(no_of_inputs)] 
plot (x)
show()
x_fft = fft.rfft(x)
plot (20*log10(abs(x_fft)))
show()
def adc(fx,fx_max,bit_num):
    output=[]
    for each_val in fx:
        dig_out_adc=''
        cmp_val = (fx_max)/2.0
        for j in arange(bit_num):
            if each_val >= cmp_val:
                dig_out_adc = dig_out_adc+'1'
                if j!=bit_num - 1:
                    cmp_val = cmp_val + (fx_max+1)*1.0/(2**(j+2))
            else:
                dig_out_adc = dig_out_adc+'0'
                if j!=bit_num - 1:
                    cmp_val = cmp_val - (fx_max+1)*1.0/(2**(j+2))
        output.append(dig_out_adc)
    return output
adc_output = adc(x,16384.0,no_of_bits)

no_of_weights = 26

def binary_to_thermometric_decoder(inputs,sources):
    value = sum([int(inputs[i])*(2**i) for i in arange(len(inputs))])
    if value == 0:
        output = 0.0
    else:
        output = sum([sources[k] for k in arange(value)])
    return output
def binary_weighted_decoder(inputs,sources):
    output = sum([int(inputs[l])*sources[len(inputs)-l-1] for l in arange(len(inputs))])
    return output

###############weights#######################
no_of_dac_bits = 14
weights = zeros(no_of_weights)
inputs_parallel = zeros(26)
for i in arange(8):
    weights[i] = 2.0**i
weights[8:11] = 256.0
weights[11:26] = 1024.0
dig_input = zeros(no_of_weights)

#print weights
dac_coefficients = [weights[b] + np.random.normal(0.0,0.1) for b in arange(no_of_weights)]
output_val = zeros(no_of_inputs)
for m in arange(no_of_inputs):
    dig_input = adc_output[m]
    output_val[m] =  binary_weighted_decoder(dig_input[0:8],dac_coefficients[0:8])
    output_val[m] += binary_to_thermometric_decoder(dig_input[8:10],dac_coefficients[8:11])
    if m%200 < 100:
        output_val[m] += binary_to_thermometric_decoder(dig_input[10:13],dac_coefficients[11:26])
    else:
        output_val[m] += binary_to_thermometric_decoder(dig_input[10:13],dac_coefficients[12:26])
        
plot (output_val)
show()
output_fft = fft.rfft(output_val,10000)
plot (20*log10(abs(output_fft)))
show()
print x_fft[50]
print output_fft[50]
