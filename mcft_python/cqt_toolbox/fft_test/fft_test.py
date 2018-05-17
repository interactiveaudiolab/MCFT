from __future__ import print_function, division

import copy
import csv

import numpy as np 
from scipy.io import loadmat
import matplotlib.pyplot as plt

'''
real_part = []
imag_part = []

with open('real_part.csv') as csvfile:
	reader = csv.reader(csvfile)
	for row in reader:
		real_part.append(float(row[0]))

with open('imag_part.csv') as csvfile:
	reader = csv.reader(csvfile)
	for row in reader:
		imag_part.append(float(row[0]))

matfft = np.zeros(len(real_part),dtype=np.complex128)
for i in range(matfft.size):
	total = real_part[i] + 1j*imag_part[i]
	matfft[i] = total

signal = loadmat('signal.mat')
signal = signal['f']

pyfft = np.fft.fft(signal,axis=0)
'''

samp_rate = 16000 # sample rate
time_vec = np.arange(0,1,1/samp_rate) # time vector
pysignal = np.cos(2*np.pi*440*time_vec)
pyfft = np.fft.fft(pysignal)


matlab = loadmat('cosine.mat')
matsignal = matlab['signal']
matfft = matlab['matfft']

sigdiff = np.sum(matsignal - pysignal)
fftdiff = np.sum(matfft - pyfft)

print(sigdiff,fftdiff)

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(abs(matfft),'b:',label='matfft')
plt.legend()
ax2.plot(abs(pyfft),'g-',label='pyfft')
plt.legend()
plt.show()