from numpy import *
from pylab import *
from scipy import signal
from numpy import fft
from random import randint
import numpy as np

no_of_bits = 8  # Number of Bits
base=2     #radix
levels = base**no_of_bits       # Total Number of Quantized Levels
Vm = 8                     # Max Voltage
alpha = 0.001                     #higher order coefficient
order = 2                       #order of non linearity
################Sample and Hold Circuit###############
dig_out_adc = zeros(no_of_bits)
no_of_inputs = 10000
#################ADC successive approximation#########
x = [1.0 for t in arange(no_of_inputs)]

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

#############Single input ADC feedback############
def adc_single_input(fx,fx_max,bit_num):
    dig_out_adc=''
    cmp_val = (fx_max)/2.0
    output=''
    for j in arange(bit_num):
        if fx >= cmp_val:
            dig_out_adc += '1'
            if j!=bit_num - 1:
                cmp_val += (fx_max+1)*1.0/(2**(j+2))
        else:
            dig_out_adc += '0'
            if j!=bit_num - 1:
                cmp_val -= (fx_max+1)*1.0/(2**(j+2))
    output+=dig_out_adc
    return output

###############Adaptive LMS filter################
no_of_weights = no_of_bits
error = zeros(no_of_inputs)
def adaptive_lms_filter(x_input,weight, desired_d):
    L = len(weight)                        
    beta = 0.0001
    #output = sum([int(x[u])*weight[u] for u in arange(L)])
    error = desired_d - x_input
    #print error,weight,x_input,beta
    for p in arange(L):
        weight[p] = weight[p] + 2*beta*error*x_input
    return (error,weight)
weights = zeros(no_of_weights)
array_mem = zeros(no_of_weights)
radix_set = [2**(no_of_bits-h-1) for h in arange(no_of_bits)]
#print weights
for m in arange(no_of_inputs):
    dig_input=adc_output[m]
    voltage_dac_input  = sum([int(dig_input[q])*radix_set[q] for q in arange(len(dig_input))])
    modified_dig_input = [int(dig_input[r])*(1.0+ weights[r]) for r in arange(len(dig_input))]
    #print modified_dig_input
    voltage_dac_output = sum([modified_dig_input[f]*radix_set[f] for f in arange(len(dig_input))])
    voltage_non_linear = voltage_dac_output + alpha*(voltage_dac_output**order)
    print voltage_dac_input,modified_dig_input,voltage_dac_output,voltage_non_linear,weights,voltage_dac_output
    temp_array = adaptive_lms_filter(voltage_non_linear/100.0, weights, voltage_dac_input/100.0)
    weights = temp_array[1]
    error[m] = temp_array[0]
    #print temp_array
    ################Adaptive Filter################
    '''
    array_mem[0] = actual_voltage
    array_mem = array(array_mem)
    array_mem[1:] = array_mem[0:len(array_mem)-1]
    '''
#array_2 = array_2/max(array_2)*(levels - 1)
#print dnlandinl(levels,array_2,d,total_steps)
#plot(array_2)
#show()
#plot (weights)
#show()
plot (error)
show()
