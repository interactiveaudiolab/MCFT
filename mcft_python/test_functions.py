import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import time as time
from scipy.fftpack import fftshift
import cqt_toolbox.cqt
from mcft_toolbox.spectro_temporal_fbank import *
from mcft_toolbox.mcft import *
from mcft_toolbox.inv_mcft import *

############################# Compute filter default centers ####################################
scale_params = (2,100,12)
rate_params = (2,200,10)
scale_ctrs,rate_ctrs = filt_default_centers(scale_params,rate_params)

############################# Compute a scale-rate filter #######################################
scale_ctr = 1
rate_ctr = 4

# scale parameters
samprate_spec = 24
scale_filt_len = 3*samprate_spec
type = 'bandpass'

scale_params = {'scale_filt_len':scale_filt_len,'samprate_spec':samprate_spec,'type':type}

# rate parameters
samprate_temp = 150
rate_filt_len = 1*samprate_temp
time_const = 3.5
type = 'bandpass'

rate_params = {'time_const':time_const,'rate_filt_len':rate_filt_len,'samprate_temp':samprate_temp,'type':type}

filt_tf_up, filt_sr_up = gen_filt_scale_rate(scale_ctr,rate_ctr,scale_params,rate_params,'up')
filt_tf_down, filt_sr_down = gen_filt_scale_rate(scale_ctr,rate_ctr,scale_params,rate_params,'down')


filt_tf_up = fftshift(filt_tf_up,axes=0)
filt_tf_down = fftshift(filt_tf_down,axes=0)

# generate frequency and time vectors for plotting
num_f,num_t = np.shape(filt_tf_up)
freq_vec = np.arange(num_f)/samprate_spec
time_vec = np.arange(num_t)/samprate_temp

# plot the filters
plt.figure()
plt.subplot(211)
plt.pcolormesh(time_vec,freq_vec,np.real(filt_tf_up))
plt.ylabel('octave number')
plt.xlabel('time (sec)')
plt.title('upward ripple')

plt.subplot(212)
plt.pcolormesh(time_vec,freq_vec,np.real(filt_tf_down))
plt.ylabel('octave number')
plt.xlabel('time (sec)')
plt.title('downward ripple')
plt.show()

# compare python and matlab results
mat_results = loadmat('mat_files/matlab_results_filt.mat')
err_up = np.sum(np.abs(filt_tf_up - mat_results['h_up']))
err_down = np.sum(np.abs(filt_tf_down - mat_results['h_down']))

#################################### Compute a scale-rate filterbank #####################################
# generate a complex signal for phase modulation
# make time vectors for constructing signals
samp_rate = 16000 # sample rate
time_vec = np.mat(np.arange(0,1,1/samp_rate)).T # time vector

# cqt parameters
fmin = 27.5*2**(0/12)
fmax = 27.5*2**(87/12)
fres = 24 # bins per octave
gamma = 10

# generate a signal and compute its cqt
signal = np.cos(2*np.pi*220*time_vec)+np.cos(2*np.pi*440*time_vec)+np.cos(2*np.pi*880*time_vec)
Xcq = cqt(signal, fres, samp_rate, fmin, fmax, gamma=gamma, window_name='hann')
comp_specgram = Xcq['cqt']

plt.pcolormesh(np.abs(comp_specgram))

# filterbank parameters
scale_ctrs = 2.**np.arange(-2,4)
rate_ctrs = 2.**np.arange(-2,7)
nfft_scale,nfft_rate = np.shape(comp_specgram)
filt_params = {'samprate_spec':fres,'samprate_temp':nfft_rate,'time_const':2}

start = time.time()
fbank_tf_domain, fbank_sr_domain = gen_fbank_scale_rate(scale_ctrs,rate_ctrs,nfft_scale,nfft_rate,filt_params,comp_specgram=comp_specgram)
end = time.time()
print('computation time:',end-start)

plt.subplot(221)
plt.pcolormesh(fftshift(np.abs(fbank_tf_domain[2,2,:]),axes=0))
plt.subplot(222)
plt.pcolormesh(fftshift(np.abs(fbank_tf_domain[2,-3,:]),axes=0))
plt.subplot(223)
plt.pcolormesh(fftshift(np.abs(fbank_sr_domain[2,2,:])))
plt.subplot(224)
plt.pcolormesh(fftshift(np.abs(fbank_sr_domain[2,-3,:])))

# compare python and matlab results
mat_results = loadmat('mat_files/matlab_results_fbank.mat')

err_spec = np.sum(np.abs(comp_specgram - mat_results['comp_specgram']))

err_tf = np.abs(fbank_tf_domain - mat_results['h_out'])
err_sr = np.abs(fbank_sr_domain - mat_results['H_out'])

err_tf = np.sum(np.abs(fbank_tf_domain - mat_results['h_out']))
err_sr = np.sum(np.abs(fbank_sr_domain - mat_results['H_out']))

############################# Test cqt-to-mcft #######################################

sig_cqt = comp_specgram
fbank_scale_rate = fbank_sr_domain

start = time.time()
mcft_out = cqt_to_mcft(sig_cqt,fbank_scale_rate)
end = time.time()
print('computation time:',end-start)

mat_results = loadmat('mat_files/cqt_to_mcft.mat')

err_mcft = np.sum(np.abs(mcft_out - mat_results['mcft_out']))

############################# Test mcft #######################################

cqt_params_in = {'samprate_sig':samp_rate,'fmin':fmin,'fmax':fmax,'fres':fres,'gamma':gamma}

start = time.time()
mcft_out, cqt_params_out, fbank_sr_domain,scale_ctrs,rate_ctrs = mcft(signal, cqt_params_in, filt_params_in=None,del_cqt_phase=0)
end = time.time()
print('computation time:',end-start)

mat_results = loadmat('mat_files/mcft.mat')

err_mcft = np.sum(np.abs(mcft_out - mat_results['mcft_out']))
err_H = np.sum(np.abs(fbank_sr_domain - mat_results['H']))

plt.pcolormesh(fftshift(np.abs(fbank_sr_domain[0,0,:,:])))

############################# Test mcft_to_cqt #######################################

start = time.time()
sig_cqt_rec = mcft_to_cqt(mcft_out,fbank_sr_domain)
end = time.time()
print('computation time:',end-start)

mat_results = loadmat('mat_files/mcft_to_cqt.mat')

err_Xhat = np.sum(np.abs(sig_cqt_rec - mat_results['X_hat']))

############################# Test inv_mcft #######################################

start = time.time()
signal_rec = inv_mcft(mcft_out, cqt_params_out, fbank_sr_domain)
end = time.time()
print('computation time:',end-start)

mat_results = loadmat('mat_files/inv_mcft.mat')

err_xhat = np.sum(np.abs(signal_rec - mat_results['x_hat']))


def norm2(x):
    return np.sqrt(np.sum(np.square(x)))

rec_err = 20*np.log10(norm2(signal-signal_rec)/norm2(signal))
