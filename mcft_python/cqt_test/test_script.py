from __future__ import print_function, division

import time
import librosa
import matplotlib.pyplot as plt
import numpy as np
from cqt import cqt
from icqt import icqt

#### TODO: Add python3 type hints to function definitions

# make time vectors for constructing signals   
samp_rate = 16000 # sample rate
time_vec = np.mat(np.arange(0,1,1/samp_rate)).T # time vector

# cqt parameters
fmin = 27.5*2**(0/12)
fmax = 27.5*2**(87/12)
fres = 24 # bins per octave
gamma = 10

# Construct 4 unique signals
signal1 = np.cos(2*np.pi*220*time_vec)+np.cos(2*np.pi*440*time_vec)+np.cos(2*np.pi*880*time_vec)
signal2, sr2 = librosa.core.load('audio_files/original/chirps.wav')
signal3, sr3 = librosa.core.load('audio_files/original/trombone.wav')

total_sig_length = (len(signal1)/samp_rate)+(len(signal2)/sr2)+(len(signal3)/sr3)

# Set start time for benchmarking
start = time.clock()

# Compute cqts of all four signals
Xcq1 = cqt(signal1, fres, samp_rate, fmin, fmax, gamma=gamma, window_name='boxcar')
Xcqt1 = Xcq1['cqt']

Xcq2 = cqt(signal2, fres, sr2, fmin, fmax, gamma=gamma, window_name='cos')
Xcqt2 = Xcq2['cqt']

Xcq3 = cqt(signal3, fres, sr3, fmin, fmax, gamma=gamma, window_name='gauss')
Xcqt3 = Xcq3['cqt']

# Check how long cqt computation took
cqt_time = time.clock()

# Take the icqt before plotting for timing purposes
re_signal1, gd = icqt(Xcq1)
re_signal2, gd = icqt(Xcq2)
re_signal3, gd = icqt(Xcq3)

# Clock the icqt time
icqt_time = time.clock() - cqt_time

# Write the reconstructed signals to file
librosa.output.write_wav('audio_files/reconstructed/re_synth.wav',signal1,samp_rate)
librosa.output.write_wav('audio_files/reconstructed/re_chirps.wav',signal2,sr2)
librosa.output.write_wav('audio_files/reconstructed/re_trombone.wav',signal3,sr3)

# Compute reconstruction error and output times
def norm2(x):
    return np.sqrt(np.sum(np.square(x)))

print('Total length (all) signals:',total_sig_length)
print('Time to compute cqt in secs:',cqt_time,'\nTime to compute icqt in secs:',icqt_time,'\n')

print('Signal 1 reconstruction error:',20*np.log10(norm2(re_signal1-signal1)/norm2(signal1)))
print('Signal 2 reconstruction error:',20*np.log10(norm2(re_signal2-signal2)/norm2(signal2)))
print('Signal 3 reconstruction error:',20*np.log10(norm2(re_signal3-signal3)/norm2(signal3)))

# Plot all 4 cqts together
plt.figure(figsize=(12,4))
plt.subplot(131)
plt.pcolormesh(abs(Xcqt1))
plt.subplot(132)
plt.pcolormesh(abs(Xcqt2))
plt.subplot(133)
plt.pcolormesh(abs(Xcqt3))
plt.show()