from __future__ import division
from cqt import cqt
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import librosa

# generate the signal    
samp_rate = 16000 # sample rate
time_vec = np.arange(0,1,1/samp_rate) # time vector
num_harmonics = 4
harmonic_freqs = np.arange(1,5)*(27.5*2**((49-1)/12))
mod_amp = np.arange(1,5)*20
mod_freq = 5

time_vec.shape = (1,len(time_vec))
harmonic_freqs.shape = (num_harmonics,1)
mod_amp.shape = (len(mod_amp),1)

# harmonic signal without frequency modulation
x1 = np.cos(2*np.pi*harmonic_freqs*time_vec)
x1 = np.sum(x1,axis=0)

# harmonic signal with frequency modulation
x2 = np.cos(2*np.pi*harmonic_freqs*time_vec+mod_amp*np.cos(2*np.pi*mod_freq*time_vec))
x2 = np.sum(x2,axis=0)
x = x1+x2


# cqt of the signal 
fmin = 27.5*2**(20/12)
fmax = 27.5*2**(87/12)
fres = 48 # bins per octave
gamma = 100

pysignal = np.cos(2*np.pi*440*time_vec)
audio,sr = librosa.core.load('chirps.wav')\

Xcq = cqt(audio, fres, samp_rate, fmin, fmax,gamma=gamma)
Xcqt = Xcq['c']

matlab = loadmat('matchirps.mat')
mcqt = matlab['Xcqt']

imag_diff = np.sum(np.imag(Xcqt) - np.imag(mcqt))
real_diff = np.sum(np.real(Xcqt) - np.real(mcqt))
print(imag_diff,real_diff)

Nf,Nt=Xcqt.shape;

freq_vec_plot=Xcq['fbas']
time_vec_plot=np.linspace(0,time_vec[:,-1],Nt);

plt.pcolormesh(abs(Xcqt))
plt.show()

'''
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.pcolormesh(time_vec_plot,freq_vec_plot,Xcqt,label='Python CQT')
ax1.set_title("Python CQT")
ax2.pcolormesh(time_vec_plot,freq_vec_plot,mcqt,label='MATLAB CQT')
ax2.set_title("MATLAB CQT")
plt.show()
'''