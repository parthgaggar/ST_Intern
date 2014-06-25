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
#ideal_step_size = Vm*1.0/(levels)   #voltage step size
value = zeros(10)               
array_2 = zeros(10)
points=zeros(10)
final_array=zeros(10)
current_sources=zeros(base**no_of_bits - 1)
no_of_iterations = 1
################Sample and Hold Circuit###############
dig_out_adc = zeros(no_of_bits)
no_of_inputs = 1000
#################ADC successive approximation#########
x = [np.random.normal(5,1) for t in arange(no_of_inputs)]

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


    
#################measuring DNL and INL#################
def dnlandinl(no_of_inputs,arrays,iteration_no,total_steps):
    dnl = [[0 for x in arange(no_of_inputs)] for x in arange(no_of_iterations)]
    inl = [[0 for x in arange(no_of_inputs)] for x in arange(no_of_iterations)]
    dnl[iteration_no][0] = 0
    inl[iteration_no][0] = 0
    for a in arange(2,levels+1):
        dnl[iteration_no][a-1] = abs(((arrays[total_steps*a - 1] - arrays[total_steps*(a-1)])*1.0/ideal_step_size) - 1)
        inl[iteration_no][a-1] = abs((arrays[total_steps*a - 1] - fx[a-1])*1.0/ideal_step_size)
    #plot(inl[d])
    #show()
    dnl_array = zeros(no_of_inputs)
    mean_inl_array = zeros(no_of_inputs)
    mean_dnl_array = zeros(no_of_inputs)

    for col in arange(no_of_inputs):
        sum_num1 = 0.0
        sum_num2 = 0.0
        for row in arange(no_of_iterations):
            sum_num1+=inl[row][col]
            sum_num2+=dnl[row][col]        
        mean_inl_array[col] = sum_num1*1.0/no_of_iterations
        mean_dnl_array[col] = sum_num2*1.0/no_of_iterations

    plot (mean_inl_array)
    plot (mean_dnl_array)
    show()
    dnl_array = [max(dnl[p]) for p in arange(no_of_iterations)]
    inl_array = [max(inl[p]) for p in arange(no_of_iterations)]

    mean_dnl = mean(dnl_array)
    mean_inl = mean(inl_array)
    return (mean_dnl,mean_inl)



##############Single input ADC feedback############
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
no_of_weights = 8
weights = zeros(no_of_weights)

error = zeros(no_of_inputs)
def adaptive_lms_filter(x,weights, desired_d):
    L = len(weights)                        
    beta = 0.001
    output = sum([x[u]*weights[u] for u in arange(L)])
    error = desired_d - output
    #print error
    weights = [weights[p]+2*beta*error*x[p] for p in arange(L)]
    return (weights,error)
weights = zeros(no_of_weights)
array_mem = zeros(no_of_weights)
transposed_matrix = [[1 for x in arange(no_of_inputs)] for x in arange(levels)]
source_matrix = array(transposed_matrix).transpose()
array_2 = zeros(no_of_inputs)
radix_set = [2**(no_of_bits-h-1) for h in arange(no_of_bits)]
#print weights
for m in arange(no_of_inputs):
    dig_input=adc_output[m]
    current_sources = source_matrix[m]
    voltage_dac_input = sum( [int(dig_input[q])*radix_set[q]) for q in arange(len(dig_input))])
    voltage_dac_output = sum([int(dig_input[f])*radix_set[f]) for f in arange(len(dig_input))])
    voltage_non_linear = voltage_dac_output + alpha*(voltage_dac_output**order)
    #print voltage_non_linear
    desired_voltage = voltage_dac_input*1.0     
    actual_voltage = voltage_non_linear*1.0
    print desired_voltage,actual_voltage
    ################Adaptive Filter################
    array_mem[0] = actual_voltage
    temp_array = adaptive_lms_filter(array_mem, weights, desired_voltage)
    weights = temp_array[0]
    error[m] = temp_array[1]
    array_mem = array(array_mem)
    array_mem[1:] = array_mem[0:len(array_mem)-1]
#array_2 = array_2/max(array_2)*(levels - 1)
#print dnlandinl(levels,array_2,d,total_steps)
plot(array_2)
show()
plot (error)
show()
