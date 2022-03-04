#!/usr/bin/python

#Routine for calculating power-spectrum of a time-series,
#given as an ASCII file with two or three columns, (time,observation,sigma).

import sys
from optparse import OptionParser
from types import *
from matplotlib import pyplot
from matplotlib.pyplot import xlabel, ylabel, savefig, ion, ioff
from math import *
from numpy import arange, loadtxt, zeros, ones, cos, sin, median, ediff1d, array, where, max as amax, empty, arctan2, arctan
from numexpr import evaluate as neeval

#DEFINING FUNCTIONS USED IN ALGORITHM

#Define function for calculating the power-spectrum
def powerspec(freqs, obs, time, weight):
    #Defining variables for calculations

    n = freqs.size

    #Making matrix with freq*time
    prods = (array([time]).repeat(n, axis=0)).T * freqs

    #Making matrix with weights
    weights = (array([weight]).repeat(n, axis=0)).T

    #Taking sine/cosine value of all freq*time
    sines, cosines = neeval("sin(prods)"), neeval("cos(prods)")

    #Creating matrix with all obs values
    obsT = array([obs]).repeat(n, axis = 0).T
    wobsT = neeval("obsT * weights")
    wsines, wcosines = neeval("sines * weights"), neeval("cosines * weights")

    #Calculating  weighted s, c, ss, cc, sc for all frequencies
    s = neeval("(sines * wobsT)").sum(axis = 0)
    c = neeval("(cosines * wobsT)").sum(axis = 0)
    ss = neeval("(sines * wsines)").sum(axis = 0)
    cc = neeval("(cosines * wcosines)").sum(axis = 0)
    sc = neeval("(sines * wcosines)").sum(axis = 0)

    #Calculating alpha and beta for all frequencies
    _preab = neeval("ss*cc-sc*sc")
    alpha = neeval("(s*cc - c*sc)/_preab")
    beta = neeval("(c*ss-s*sc)/_preab")

    return alpha, beta

def chunkpowexp(freq, data, time, weight):
    if freq.size < 1000:
       result = array(powerspec(afreq, data, time, weight))
    else:
       result = empty([2,len(freq)])
       for i in range(0, freq.size, 1000):
            result[0, i:i+1000], result[1, i:i+1000] = powerspec(afreq[i:i+1000], data, time, weight)

    return result[0], result[1], (result * result).sum(axis = 0)

#Making filter function
def filter(index,freq,time,params):
    #Making index values
    m = time.size

    #Making matrices needed for calculation
    alpha = params[0,index].repeat(m, axis=0).T
    beta = params[1,index].repeat(m, axis=0).T
    prods = array([freq[index]]).repeat(m, axis=0).T * time
    #Taking sine/cosine for all frequencies in range
    sines, cosines = neeval("sin(prods)"), neeval("cos(prods)")

    #Making matrix for all values in final sum
    summatrix = neeval("alpha * sines + beta * cosines")
    #Calculating sum for all times in timeserie
    filtersum = summatrix.sum(axis = 0)
    return filtersum


#-------------------------MAIN PROGRAM STARTS HERE--------------------------------

#General usage of program:
usage = "usage: %prog [options] filepath"
parser = OptionParser(usage=usage)
parser.add_option("-d","--delim", dest="dl", default = None, help = 'Sets the file delimiter, default is any whitespace.')
parser.add_option("-s","--skip", type = "int", dest='n', default = 0, help = 'Skip first N lines of file, default is 0.')
parser.add_option("-t","--time",type = "string", dest='t', default = 's', help = 'Time measured in days (d), hours (h), minutes (m) or seconds (s). Default is seconds')
(options,args)=parser.parse_args()

#Check if filepath is given
if len(args) < 1:
   print "Error: No file-path specified. Type -h for more information."
   sys.exit(0)

path = args[0]

#Sorts data into new variables
data = loadtxt(path,delimiter=options.dl, skiprows=options.n)
time = data.T[0]
if options.t == 'd':
        #Changing from days to seconds
        time = time * 24 * 60 * 60
if options.t == 'h':
        #Changing from hours to seconds
        time = time * 60 * 60
if options.t == 'm' :
        #Changing from minutes to seconds
        time == time * 60
obs = data.T[1]


#Loading weights from data file
#Set equal to one if no weights present
if data.shape[1] < 3:
        weight = ones(len(time))
else :
        weight = data.T[2]**(-2.)

#weight = ones(len(time))

#Making harmonic oscillator for testing of routine
#time = arange(0,30,.01)
#obs =sin(2.*pi*1.5*time)
#weight = ones(len(time))
#obs = 0.7 * sin(2. * pi * 1.5 * time)+sin(2.*pi*1.8*time+2*pi) + 0.15 * sin(2*pi*0.6 * time+pi/4) + 0.9 * sin(2*pi*2.1*time-0.23*2*pi) + 0.5 * sin(2*pi*1.3*time+0.12*2*pi)

#Making simulation array for selected range of solar oscillations
#obs = 0.32994*sin(2. * pi * 0.0031587 * time) + 0.33046*sin(2.*pi*0.00316817*time) + 0.38569*sin(2*pi*0.00323285*time)

#Calculate minimum frequency and resolution:
deltat = time[len(time)-1]-time[0]
nu_min = 1./deltat
delta_nu = 1./(10.*deltat)

#Calculating Nyquist frequency of time-series:
nu_ny = 1./(2*median(ediff1d(time)))
print 'Please enter frequency range and frequency resolution for the power-spectrum.'
print 'Suggestion: Freq. range: [',nu_min, nu_ny,'], res:', delta_nu

#Promt user for frequency interval and resolution
fstart, fend, res =  input('Please enter Freq. start, Freq end, Freq. resolution: ')
freq = arange(float(fstart),float(fend),res)

#Subbstracting mean-value from timeseries
obsm = obs - obs.mean()

# We want angular frequency
afreq = freq * 2 * pi

#Calculating powerspectrum:
#Splits data in chunks if necessary
alpha, beta, power = chunkpowexp(afreq, obsm, time, weight)

#Showing powerspectrum
pyplot.ion()
pyplot.figure()
pyplot.plot(freq,power)
pyplot.xlabel('Frequency')
pyplot.ylabel('Power')
pyplot.title('Powerspectrum')

#Calculating window-function for observations
wfreq1 = 2 * pi * input('Please enter window frequency: ') * time
windowp = 0.5 * chunkpowexp(afreq, sin(wfreq1) + cos(wfreq1),
                              time, weight)[-1]

#Implementing high-pass filter
checkh = raw_input('Do high-pass filtering (y/n)? ') == 'y'
if checkh:
    hfreq = float(raw_input('High-pass frequency: '))
    hindex = where(freq<hfreq)
    index = where(freq>hfreq)
    highpass = obs - filter(hindex,afreq,time,param)/windowp.sum()
    highpassm = highpass - highpass.mean()
    fpower = chunkpowexp(afreq, highpassm, time, weight)[-1]

#Implementing low-pass filter
checkl = raw_input('Do low-pass filtering (y/n)? ') == 'y'
if checkl:
    lfreq = float(raw_input('Low-pass frequency: '))
    index = where(freq<lfreq)
    lowpass = filter(index, afreq, time, param)/windowp.sum()
    lpassm = lowpass - lowpass.mean()
    fpower = chunkpowexp(afreq, lpassm, time, weight)[-1]

#Implementing band-pass filter
checkb = raw_input('Do band-pass filtering (y/n)? ') == 'y'
if checkb:
    bstart, bend = input('Please enter band-pass range. Freq. start, Freq. end: ') 
    index = where(freq == (freq-(freq<bstart)-(freq>bend)))
    bandpass = filter(index, afreq, time, param)/windowp.sum() 
    bpassm = bandpass - bandpass.mean()
    fpower = chunkpowexp(afreq, bpassm, time, weight)[-1]

#Showing filtered spectrum
if checkh or checkl or checkb:
    pyplot.close()
    pyplot.plot(freq[index],fpower[index])
    pyplot.xlabel('Frequency')
    pyplot.ylabel('Power')
    pyplot.title('Filtered spectrum')

#Setting observed spectrum to filtered spectrum
if checkh:
    obsm=highpassm
elif checkl:
    obsm=lpassm
elif checkb:
    obsm=bpassm

#Implementing CLEAN algorithm.

checkclean = raw_input('CLEAN spectrum? (y/n) ') == 'y'
if checkclean:
    #Prompting user for CLEANing input
    cln = int(raw_input('Number of CLEAN iterations: ')) 
    clstart, clend, cres = input('Please enter CLEAN Freq. start, Freq. end, resolution: ')
    clfreq = arange(clstart, clend, cres)
    aclfreq = clfreq * 2 * pi
    robs = obsm
    cleanparam = empty([cln,3])
    #Calculating powerspectrum in CLEAN range
    alpha, beta, cpower = chunkpowexp(aclfreq, robs, time, weight)
    #CLEANing spectrum, lower resolution.
    for i in range(cln):
        m = where(cpower == amax(cpower))
        cleanparam[i,0], cleanparam[i,1],cleanparam[i,2] =  clfreq[m][0], cpower[m][0],arctan(cparam[1,m][0]/cparam[0,m][0])
        _aclfreq = aclfreq[m] * time
        robs = robs - (alpha[m] * sin(_aclfreq)) - (beta[m] * cos(_aclfreq))
        cpower = chunkpowexp(aclfreq, robs, time, weight)[-1]
    #Recalculating CLEAN spectrum in full resolution
    cpower = chunkpowexp(afreq, robs, time, weight)[-1]
else :
    cpower = power

#Setting path and filename for results-file

#path = '/home/aot/time-series-analysis/keplerwindow.txt'
#path = '/home/aot/time-series-analysis/phasetest.txt'
path1=path.split('.')


respath = path1[0] + '-power.dat'
respath2 = path1[0] + '-clean.dat'
respath3 = path1[0] + '-results.dat' 
respath4 = path1[0] + '-filter.dat'
respath5 = path1[0] + '-window.dat'

if checkh:
    fpower=hpower
elif checkl:
    fpower=lpower
elif checkb:
    fpower=bpower

if checkh and checkl and checkb:
    fpower = power
    index = where(freq == freq)

pyplot.close()
pyplot.ioff()
pyplot.figure(1)
pyplot.subplot(131)
pyplot.plot(freq[index],fpower[index])
pyplot.xlabel('Frequency')
pyplot.ylabel('Power')
pyplot.title('Powerspectrum')
pyplot.subplot(132)
pyplot.plot(freq[index],cpower[index])
pyplot.xlabel('Frequency')
pyplot.ylabel('Power')
pyplot.title('CLEANed spectrum')
pyplot.subplot(133)
pyplot.plot(freq,windowp)
pyplot.xlabel('Frequency')
pyplot.ylabel('Power')
pyplot.title('Windowfunction')
pyplot.show()

def savedata(path, data, descr):
    with open(path, 'w') as result:
        for (f, p) in zip(freq, data):
            print >> result, "%g\t%g" % (f, p)
        print ('%s saved as: %s' % (descr, respath))

#Writing data into new files
savedata(respath, power, 'Power spectrum')
savedata(respath5, windowp, 'Window function')
if checkclean:
    savedata(respath2, cpower, 'CLEANed spectrum')
    with open(respath3, 'w') as fresults:
        for (f,p,ph) in cleanparam:
            print >> fresults, "%g\t%g\t%g" % (f,p,ph)
        print 'Extracted frequencies saved as: ', respath3

if checkh or checkl or checkb:
    savedata(respath4, fpower, 'Filtered powerspectrum')
