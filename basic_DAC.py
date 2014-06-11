from numpy import *
from pylab import *
from numpy import fft
from random import randint

no_of_bits = 3  # Number of Bits 
base=2     #radix
levels = base**no_of_bits  # Total Number of Quantized Levels
Vm = levels             # Max Voltage
time_init = 0       #start time
time_final = time_init + 1     #end time
time_step = 0.1    #time steps
rvalue = 1          #load resistor
cvalue = 0.001        #load capacitance
count = 0
#dig_input=arange(10)
value = zeros(10)
array_2 = zeros(10)
points=zeros(10)
final_array=zeros(10)
fig=figure()
################Sample and Hold Circuit###############
dig_out_adc = zeros(no_of_bits)
sample_length = 1.0
fx_max = 7.0
sample_duration = 0.001
#samples = sample_length/sample_duration
t = np.arange(0 ,sample_length, sample_duration)
f=10.0
fx = fx_max/2.0+fx_max/2.0 * sin(2*pi*f*t/sample_length)
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
        cmp_val = (fx_max+1)/2.0
        for j in arange(bit_num):
            #print cmp_val
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

adc_output = adc(fx,fx_max,no_of_bits)
#adc_output = ['000','001','010','011','100','101','110','111']
no_of_inputs = len(adc_output) #input count
                    #######DAC starts here########

dig_input=adc_output
def current_cell(x,radix):
    length = len(x)
    output = sum([int(x[f])* (radix**(length-f-1))for f in arange(length)])
    return output
    
def first_order_out_rc(step_size, time_init, time_final, time_step, rvalue, cvalue):
    sample_points = arange(time_init, time_final, time_step)
    val= step_size*(1 - exp (-1/(rvalue*cvalue)*sample_points))
    return val

for count in arange(no_of_inputs):
    voltage = current_cell(dig_input[count],base)  #output voltage
    normalized_voltage = 1.0*voltage/(base**no_of_bits)  #normalized voltage
    desired_voltage = normalized_voltage * Vm     #Final Output voltage
    points = arange(time_init, time_final, time_step)
    arr=array_2
    tmpval= arr[len(arr)-1]
    array_1 = first_order_out_rc(desired_voltage-tmpval, time_init+1, time_final+1, time_step, rvalue, cvalue)
    array_2 = [tmpval + array_element for array_element in array_1]
    if count!=0:
        array_2= concatenate((arr,array_2), axis=0)

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
