from __future__ import division
from cqt import cqt
import numpy as np
import matplotlib.pyplot as plt

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


Xcq = cqt(x, fres, samp_rate, fmin, fmax,gamma=gamma)
Xcqt = Xcq['c']

import pdb; pdb.set_trace()

Nf,Nt = Xcqt.shape
#fvec = Xcq['fbas']
#tvec = np.linspace(0,time_vec[-1],Nt)

plt.imshow(Xcqt.astype(np.float64), cmap='hot')
plt.show()