from __future__ import division

import librosa
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from cqt import cqt
from icqt import icqt

# generate the signal    
samp_rate = 16000 # sample rate
time_vec = np.arange(0,1,1/samp_rate) # time vector

time_vec.shape = (1,len(time_vec))


# cqt of the signal 
fmin = 27.5*2**(20/12)
fmax = 27.5*2**(87/12)
fres = 48 # bins per octave
gamma = 100

pysignal = np.cos(2*np.pi*440*time_vec)
audio,sr = librosa.core.load('chirps.wav')

Xcq = cqt(audio, fres, samp_rate, fmin, fmax,gamma=gamma)
Xcqt = Xcq['c']

re_audio, fbank = icqt(Xcq)

# # # # # 
# # This is used for visualizing the CQT
# imag_diff = np.sum(np.imag(Xcqt) - np.imag(mcqt))
# real_diff = np.sum(np.real(Xcqt) - np.real(mcqt))
# print(imag_diff,real_diff)

# Nf,Nt=Xcqt.shape;

# freq_vec_plot=Xcq['fbas']
# time_vec_plot=np.linspace(0,time_vec[:,-1],Nt);

# plt.pcolormesh(abs(Xcqt))
# plt.show()
# # # # #