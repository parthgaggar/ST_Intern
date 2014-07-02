from numpy import *
from pylab import *
from scipy import signal
from numpy import fft
from random import randint
import numpy as np

no_of_bits = 1  # Number of Bits
base=2     #radix
Vm = 8                          # Max Voltage
alpha = 0.001                   #higher order coefficient
order = 2                       #order of non linearity
################Sample and Hold Circuit###############
dig_out_adc = zeros(no_of_bits)
no_of_inputs = 16000
#################ADC successive approximation#########
x = [np.random.normal(0,1) for t in arange(no_of_inputs)] #signal with mean = 0 and std.dev = 1
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
adc_output = adc(x,max(x),no_of_bits)

no_of_weights = 26 
error = zeros(no_of_inputs)

################LMS#########################
def adaptive_lms(x_input,weight, desired_d,beta):
    L = len(weight)                        
    y_output = sum([x_input[a]*weight[a] for a in arange(no_of_weights)])
    error = (desired_d - y_output)/1.0
    #print error,weight,x_input,beta
    for p in arange(L):
        weight[p] += 2*beta*error*x_input[p]
    return (error,weight)
###############weights#######################
no_of_dac_bits = 14
weights = zeros(no_of_weights)

for i in arange(8):
    weights[i] = 2.0**i
weights[8:11] = 256.0
weights[11:26] = 1024.0
dig_input = zeros(no_of_weights)
#print weights


dac_coefficients = [weights[b] + np.random.normal(0.0,0.01) for b in arange(no_of_weights)]

for m in arange(no_of_inputs):
    dig_input[0] = adc_output[m]
    dac_out = sum([dig_input[q]*dac_coefficients[q] for q in arange(no_of_weights)])
    tmp_array = adaptive_lms(dig_input,weights,dac_out,0.01)
    error[m] = tmp_array[0]
    weights = tmp_array[1]
    dig_input[1:] = dig_input[0:no_of_weights-1]
    dig_input[0] = 0.0

plot (error)
show()
