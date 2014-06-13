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
mu = 1                          #mean of radix
sigma = 0.0                     #std dev of radix
alpha = 0.001                     #higher order coefficient
order = 2                       #order of non linearity
ideal_step_size = Vm*1.0/(levels)   #voltage step size
value = zeros(10)               
array_2 = zeros(10)
points=zeros(10)
final_array=zeros(10)
current_sources=zeros(base**no_of_bits - 1)
no_of_iterations = 1
################Sample and Hold Circuit###############
dig_out_adc = zeros(no_of_bits)
sample_length = 1
fx_max = Vm - 1
sample_duration = 0.001
t = np.arange(0 ,sample_length, sample_duration)
f=10.0
fx = fx_max/2.0+fx_max/2.0 * sin(2*pi*f*t/sample_length)
#fx = t
#plot(fx)
#show()

#Input SINAD
Fkk = fft.rfft(fx)
SignalPow=(abs(Fkk[0])**2 + abs(Fkk[f])**2)
PowerSq=abs(Fkk)**2
PowerSq[0] = 0
PowerSq[f]=0
NoisenSigPower=sum(PowerSq)
sinad_input = 10*log10(SignalPow/NoisenSigPower)
print "SINAD_input:",sinad_input

#Input THD
m=arange(0,len(PowerSq))
distortion_input = [PowerSq[i] for i in m if i%f==0]
harmonic_distortion_input = sum(distortion_input)
thd_input = 10*log10(harmonic_distortion_input/SignalPow)
print "THD_input:",thd_input

#Input SNR
noise_input = [PowerSq[i] for i in m if i%f!=0]
noise_input = sum(distortion_input)
snr_input = 10*log10(SignalPow/noise_input)
print "SNR_input:",snr_input

#Input SFDR
worst_spur = max(PowerSq)
sfdr_input = 10*log10(SignalPow/worst_spur)
print "SFDR_input",sfdr_input

#Plot the signal
#plot (20*log10(abs(Fkk)))
#show()

#############ADC using successive approximation#############
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

adc_output = adc(fx, Vm-1, no_of_bits)
no_of_inputs = len(adc_output) #input count

                    #######DAC constituents########

time_init = 0                                       #start time
time_final = time_init + 1*sample_duration          #end time
time_step = 0.01*sample_duration                     #time steps
total_steps = (time_final - time_init)/time_step
rvalue = 10                                          #load resistor
cvalue = 0.001*sample_duration                      #load capacitance

def mismatch_power_calculator(current_source, radix, exponent, num_bits):
    val = 0
    skip_array =[1.0*2**x for x in arange(exponent + 1, num_bits)]
    skip_length = sum(skip_array)
    
    for i in arange(radix**exponent):
        val = val + current_source[skip_length + i]
    return val
        
    

def current_cell(current_source,x,radix,alpha,order):
    f=0
    output=0
    length = len(x)
    output = sum([int(x[f])* mismatch_power_calculator(current_source, radix, length-1-f, no_of_bits) for f in arange(length)])
    output +=alpha*(output**order)
    return output
    
def first_order_out_rc(step_size, time_init, time_final, time_step, rvalue, cvalue):
    sample_points = linspace(time_init, time_final, total_steps)
    val= step_size*(1 - exp (-1/(rvalue*cvalue)*sample_points))
    return val

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


###############Adaptive LMS filter################
no_of_weights = 12
weights = [0 for x in arange(no_of_weights)]

def adaptive_lms_filter(x,weights, desired_d):
    L = len(weights)                        #No of weights
    beta = 0.001
    output = sum([x[p]*weights[p] for p in arange(L)])
    error = desired_d - output
    for p in arange(L):
        weights[p] = weights[p]+2*beta*error*x[p]
    return (weights,error)
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
######################ADC feedback###############
no_of_bits_adc = 16
full_scale_voltage_adc = 8
dac_sampling_frequency = 1.0/sample_duration
adc_sampling_frequency = dac_sampling_frequency/10.0
adc_sampling_duration = 1.0/adc_sampling_frequency
dac_nyquist_frequency = dac_sampling_frequency/2.0
##################Butterworth filter#############
filter_order = 1
cutoff_freq_fraction = adc_sampling_frequency/dac_nyquist_frequency
b, a = signal.butter(filter_order,0.04,'low',analog = True)
w, h = signal.freqs(b, a)
plt.plot(20*log10(abs(h)))
array_mem = zeros(no_of_weights)
transposed_matrix = [[1 for x in arange(no_of_inputs)] for x in arange(levels)]
source_matrix = array(transposed_matrix).transpose()
array_2 = zeros(no_of_inputs)
array_2_filtered = zeros(no_of_inputs)

for m in arange(no_of_inputs):
    dig_input=adc_output[m]
    dec_input = sum([int(dig_input)*(2**(no_of_bits-z-1)) for z in arange(no_of_bits)])
    dec_input = dec_input/(2**no_of_bits)*(Vm)
    current_sources = source_matrix[m]
    tmpval = 0
    volt = 0
    desired_voltage = 0
    
    inl_array = zeros(no_of_iterations)
    #convert from binary to decimal
    analog_output=0
    
    voltage = current_cell(current_sources,dig_input,base,alpha,order)  #output voltage
    desired_voltage = voltage/levels * Vm*1.0     #Final Output voltage
    arr=array_2
    arr_filtered = array_2_filtered
    #plot (arr_filtered)
    #show()
    tmpval= array_2[len(array_2)-1]
    tmpval_filtered = array_2_filtered[len(array_2_filtered)-1]
    #create array_2
    array_1 = first_order_out_rc(desired_voltage-tmpval, time_init, time_final, time_step, rvalue, cvalue)
    array_2 = [tmpval + array_element for array_element in array_1]
    #plot (array_2)
    #show()
    array_2_filtered = signal.filtfilt(b, a, array_2)
    #plot (array_2_filtered)
    show()
    if m!=0:
        array_2= concatenate((arr,array_2), axis=0)
        array_2_filtered =  concatenate((arr_filtered,array_2_filtered), axis=0)
    #plot (array_2_filtered)
    #show()
    if m%(adc_sampling_duration/sample_duration) == 0:
        array_2_filtered_sampled = array_2_filtered[total_steps-1]
        #Fk_filtered = fft.rfft(array_2_filtered_sampled)
        #plot (20*log10(abs(Fk_filtered)))
        #show()
        adc_feedback_output = adc_single_input(array_2_filtered_sampled, full_scale_voltage_adc, no_of_bits_adc)
        #################Adaptive Filter################
        bin_to_dec = sum([int(adc_feedback_output[q])*(2**(no_of_bits_adc-q-1)) for q in arange(no_of_bits_adc)])
        bin_to_dec = bin_to_dec/(2**no_of_bits_adc)*(full_scale_voltage_adc)
        array_mem[0] = bin_to_dec
        temp_array = adaptive_lms_filter(array_mem, weights, dec_input)
        weights = temp_array[0]
        error = temp_array[1]
        
        array_mem = array(array_mem)
        array_mem[1:] = array_mem[0:len(array_mem)-1]
#array_2 = array_2/max(array_2)*(levels - 1)
#print dnlandinl(levels,array_2,d,total_steps)
plot(array_2)
plot (array_2_filtered)
show()
Fk_filtered = fft.rfft(array_2)
#plot (20*log10(abs(Fk_filtered)))
#show()
FF_filtered = fft.rfft(array_2_filtered)
#plot (20*log10(abs(array_2_filtered)))
#show()
#Output SINAD
Fk = fft.rfft(array_2)
m=arange(0,len(PowerSq))
SignalPower=(abs(Fk[0])**2 + abs(Fk[f])**2)
PowerSquare=abs(Fk)**2
PowerSquare[0] = 0
PowerSquare[f] = 0
noise_and_distortion = [PowerSquare[i] for i in m]
NoisenSignalPower=sum(noise_and_distortion)
sinad_output = 10*log10(SignalPower/NoisenSignalPower)
print "SINAD_output:",sinad_output

#Output THD
distortion = [PowerSquare[i] for i in m if i%f==0]
harmonic_distortion_output = sum(distortion)
thd_output = 10*log10(harmonic_distortion_output/SignalPower)
print "THD_output:",thd_output

#Output SNR
noise = [PowerSquare[i] for i in m if i%f!=0]
noise_output = sum(noise)
snr_output = 10*log10(SignalPower/noise_output)
print "SNR_output",snr_output

#Output SFDR
worst_spur = max(noise_and_distortion)
sfdr_output = 10*log10(SignalPower/worst_spur)
print "SFDR_output",sfdr_output
    
#Plot the signal
#plot(20*log10(abs(Fk)))
#show()
