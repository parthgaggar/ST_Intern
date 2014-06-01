from numpy import *
from pylab import *
from scipy import signal
from numpy import fft
from random import randint
import numpy as np

no_of_bits = 8  # Number of Bits
base=2     #radix
levels = base**no_of_bits       # Total Number of Quantized Levels
Vm = levels                     # Max Voltage
no_of_iterations = 1            #No of times the DAC is run
mu = 1                          #mean of radix
sigma = 0.0                     #std dev of radix
alpha = 0.0                   #higher order coefficient
order = 2                       #order of non linearity
transposed_matrix = [[0 for x in arange(no_of_iterations)] for x in arange(levels)]
for p in arange(levels):
    transposed_matrix[p] = [1 for i in arange(no_of_iterations)]

source_matrix = array(transposed_matrix).transpose()
ideal_step_size = Vm*1.0/(levels)   #voltage step size
value = zeros(10)               
array_2 = zeros(10)
points=zeros(10)
final_array=zeros(10)
current_sources=zeros(base**no_of_bits - 1)

################Sample and Hold Circuit###############
dig_out_adc = zeros(no_of_bits)
sample_length = 1
fx_max = Vm - 1
sample_duration = 0.001
t = np.arange(0 ,sample_length, sample_duration)
f=10.0
fx = fx_max/2.0+fx_max/2.0 * sin(2*pi*f*t/sample_length)
#fx = t
plot(fx)
show()

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
plot (20*log10(abs(Fkk)))
show()

#############ADC using successive approximation#############
def adc(fx,fx_max,bit_num):
    output=[]
    for each_val in fx:
        j=0
        dig_out_adc=''
        cmp_val = (fx_max)/2.0
        while j < bit_num:

            if each_val >= cmp_val:
                dig_out_adc = dig_out_adc+'1'
                if j!=bit_num - 1:
                    cmp_val = cmp_val + (fx_max+1)*1.0/(2**(j+2))
            
            else:
                dig_out_adc = dig_out_adc+'0'
                if j!=bit_num - 1:
                    cmp_val = cmp_val - (fx_max+1)*1.0/(2**(j+2))
            j=j+1
        output.append(dig_out_adc)
    return output

adc_output = adc(fx, fx_max, no_of_bits)
no_of_inputs = len(adc_output) #input count

                    #######DAC starts here########

    
dig_input=adc_output
time_init = 0                                   #start time
time_final = time_init + 1*sample_duration      #end time
time_step = 0.1*sample_duration                 #time steps
total_steps = (time_final - time_init)/time_step
rvalue = 1                                      #load resistor
cvalue = 0.001*sample_duration                                  #load capacitance
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
    
    for m in x:
        output += int(x[f])* mismatch_power_calculator(current_source, radix, length-1-f, no_of_bits)
        f+=1
    output +=alpha*(output**order)
    return output
    
def first_order_out_rc(step_size, time_init, time_final, time_step, rvalue, cvalue):
    sample_points = arange(time_init, time_final, time_step)
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




for d in arange(no_of_iterations):
    for i in arange(0,base**no_of_bits):
        current_sources = source_matrix[d]
    tmpval = 0
    count = 0
    volt = 0
    desired_voltage = 0
    array_2 = zeros(no_of_inputs)
    inl_array = zeros(no_of_iterations)
    for count in arange(no_of_inputs) :
    
        #convert from binary to decimal
        analog_output=0
    
        voltage = current_cell(current_sources,dig_input[count],base,alpha,order)  #output voltage
        
        normalized_voltage = 1.0*voltage/(levels)  #normalized voltage
    
        desired_voltage = normalized_voltage * Vm     #Final Output voltage
        arr=array_2
        tmpval= array_2[len(array_2)-1]
        if dig_input[count] == dig_input[count-1]:
            array_1 = first_order_out_rc(desired_voltage-tmpval, time_init+1*sample_duration, time_init+2*sample_duration, time_step, rvalue, cvalue)
        else:
            array_1 = first_order_out_rc(desired_voltage-tmpval, time_init, time_final, time_step, rvalue, cvalue)
    
        array_2 = [tmpval + array_element for array_element in array_1]
        if count!=0:
            array_2= concatenate((arr,array_2), axis=0)
    #array_2 = array_2/max(array_2)*(levels - 1)
    #print dnlandinl(levels,array_2,d,total_steps)
    plot(array_2)
    show()
    
    #Output SINAD
    Fk = fft.rfft(array_2)
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
    plot(20*log10(abs(Fk)))
    show()

######################ADC feedback###############

adc_sampling_duration = 10*sample_duration
adc_sampling_frequency = 1.0/adc_sampling_duration
filter_order = 4
cutoff_freq = adc_sampling_frequency/2.0 
b, a = signal.butter(filter_order,cutoff_freq,'low',analog = True)
w, h = signal.freqs(b, a)
plot (20*log10(abs(h)))
show()
dac_analog_output = [array_2[k-1] for k in arange(1*total_steps,len(array_2)+1,adc_sampling_duration)]
adc_feedback_output = adc(dac_analog_output, max(dac_analog_output), no_of_bits)
bin_to_dec = zeros(len(adc_feedback_output))

