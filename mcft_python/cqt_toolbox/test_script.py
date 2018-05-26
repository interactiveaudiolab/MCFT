from __future__ import division

import librosa
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

from cqt import cqt
from icqt import icqt

#### TODO: Add python3 type hints to function definitions

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

# # # # #
# This is for computing the icqt and reconstructing the audio
re_audio, fbank = icqt(Xcq)
librosa.output.write_wav('reconstruction.wav',re_audio,sr)
re_error = np.sum(abs(audio-re_audio))
# print(re_error)
# # # # #

# # # # # 
# This is used for visualizing the CQT
matlab = loadmat('matcqt')
mcqt = matlab['Xcqt']

imag_diff = np.sum(abs(np.imag(Xcqt) - np.imag(mcqt)))
real_diff = np.sum(abs(np.real(Xcqt) - np.real(mcqt)))
print(imag_diff,' absolute difference in imaginary components')
print(real_diff,' absolute difference in real components')

Nf,Nt=Xcqt.shape;

freq_vec_plot=Xcq['fbas']
time_vec_plot=np.linspace(0,time_vec[:,-1],Nt);

plt.pcolormesh(abs(Xcqt))
plt.show()
# # # # #