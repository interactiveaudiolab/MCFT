from __future__ import print_function, division

import time

import librosa
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from cqt import cqt
from icqt import icqt

#### TODO: Add python3 type hints to function definitions

# make time vectors for constructing signals   
samp_rate = 16000 # sample rate
time_vec = np.arange(0,1,1/samp_rate) # time vector

# cqt params
fmin = 27.5*2**(20/12)
fmax = 27.5*2**(87/12)
fres = 48 # bins per octave
gamma = 100

# Construct 4 unique signals
signal1 = np.cos(2*np.pi*440*time_vec)
signal2 = np.cos(2*np.pi*261.626*time_vec) + np.cos(2*np.pi*391.995*time_vec) 
signal3, sr3 = librosa.core.load('chirps.wav')
signal4, sr4 = librosa.core.load('trombone.wav')


# signal3, sr3 = librosa.core.load('40sec.wav')
# size = int(signal3.size)
# signal1 = signal3[:int(size/4)]
# signal2 = signal3[:int(size/2)]
# sr1 = sr3
# sr2 = sr3

# Set start time for benchmarking
start = time.clock()

# Compute cqts of all four signals
Xcq1 = cqt(signal1, fres, samp_rate, fmin, fmax, gamma=gamma, window_name='boxcar')
Xcqt1 = Xcq1['cqt']
#time1 = time.clock()

Xcq2 = cqt(signal2, fres, samp_rate, fmin, fmax, gamma=gamma, window_name='cos')
Xcqt2 = Xcq2['cqt']
#time2 = time.clock() - time1

Xcq3 = cqt(signal3, fres, sr3, fmin, fmax, gamma=gamma, window_name='gauss')
Xcqt3 = Xcq3['cqt']
#time3 = time.clock() - time2 - time1
#print(time1,time2,time3)

Xcq4 = cqt(signal4, fres, sr4, fmin, fmax, gamma=gamma)
Xcqt4 = Xcq4['cqt']

# Check how long cqt computation took
cqt_time = time.clock()

# Take the icqt before plotting for timing purposes
re_signal1, gd = icqt(Xcq1)
re_signal2, gd = icqt(Xcq2)
re_signal3, gd = icqt(Xcq3)
re_signal4, gd = icqt(Xcq4)

# Clock the icqt time
icqt_time = time.clock() - cqt_time

# Write the reconstructed signals to file
librosa.output.write_wav('A440.wav',signal1,samp_rate)
librosa.output.write_wav('fifth.wav',signal2,samp_rate)
librosa.output.write_wav('re_chirps.wav',signal3,sr3)
librosa.output.write_wav('re_trom.wav',signal4,sr4)

# Compute reconstruction error and output times
print('Time to compute cqt in secs:',cqt_time,'\nTime to compute icqt in secs:',icqt_time,'\n')
print('Signal 1 reconstruction error:',np.sum(abs(re_signal1-signal1)))
print('Signal 2 reconstruction error:',np.sum(abs(re_signal2-signal2)))
print('Signal 3 reconstruction error:',np.sum(abs(re_signal3-signal3)))
print('Signal 4 reconstruction error:',np.sum(abs(re_signal4-signal4)))

# Plot all 4 cqts together
fig = plt.figure(figsize=(12,4))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax1.pcolormesh(abs(Xcqt1))
ax2.pcolormesh(abs(Xcqt2))
ax3.pcolormesh(abs(Xcqt3))
ax4.pcolormesh(abs(Xcqt4))

plt.show()