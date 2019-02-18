import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from scipy.fftpack import fft2,ifft2,fftn,ifftn, fft, ifft, helper
import matplotlib.pyplot as plt
import cqt_toolbox.cqt
from mcft_toolbox.spectro_temporal_fbank import *
from mcft_toolbox.mcft import *
from mcft_toolbox.inv_mcft import *

import time as time
import cProfile

################################## example signal ##################################

fs = 8000
sig_dur = 2
t = np.arange(sig_dur*fs)/fs
x = np.cos(2*np.pi*440*t + 10*np.cos(2*np.pi*4*t))

################################## CQT ##################################

fmin = 27.5 * 2**(0/12)
fmax = 27.5 * 2**(87/12)
fres = 24
gamma = 0

cqt_pack = cqt(x, fres, fs, fmin, fmax,gamma=gamma)
mag_cqt = np.abs(cqt_pack['cqt'])

# plt.pcolormesh(np.abs(cqt))
# plt.show()

################################## Filterbank ##################################

# parameters
scale_res = 1
rate_res = 8

num_freq_bin,num_time_frame = np.shape(mag_cqt)
nfft_scale, nfft_rate = num_freq_bin, num_time_frame

samprate_spec = fres
samprate_temp = np.floor(num_time_frame / sig_dur)

scale_params = (scale_res,nfft_scale,samprate_spec)
rate_params = (rate_res,nfft_rate,samprate_temp)

time_const = 1
filt_params = {'samprate_spec':samprate_spec, 'samprate_temp':samprate_temp, 'time_const':time_const}

# default centers
scale_ctrs, rate_ctrs = filt_default_centers(scale_params,rate_params)

print('scale center: ', scale_ctrs, '\n\nrate centers', rate_ctrs)
print(len(scale_ctrs), len(rate_ctrs))

# filterbank
start = time.time()
fbank_sr_domain = gen_fbank_scale_rate(scale_ctrs,rate_ctrs,nfft_scale,nfft_rate,filt_params,
                                                        fbank_out_domain='sr')
stop = time.time()

print('computation time: ',stop-start)
print('fbanck dimensions: ',fbank_sr_domain.shape)

################################## MCFT ##################################

start = time.time()
# cp = cProfile.Profile()
# cp.enable()

mcft_out = cqt_to_mcft(mag_cqt,fbank_sr_domain,mcft_out_domain='tf')

# cp.disable()
# cp.print_stats()
stop = time.time()
print('computation time: ',stop-start)
print('mcft dimensions: ', mcft_out.shape)


################################## MCFT ##################################

# start = time.time()
# # cp = cProfile.Profile()
# # cp.enable()
#
# mcft_out_vec = cqt_to_mcft_vec(mag_cqt,fbank_sr_domain,mcft_out_domain='tf')
#
# # cp.disable()
# # cp.print_stats()
# stop = time.time()
# print('vec computation time: ',stop-start)
# print('mcft dimensions: ', mcft_out_vec.shape)
#
# print('err:', np.sum(np.abs(mcft_out-mcft_out_vec)))

# mcft_mean = np.mean(np.mean(np.abs(mcft_out),3),2)
# mcft2_mean = np.mean(np.mean(np.abs(mcft_out)**2,3),2)
#
# mcft_mean_vec = np.mean(np.mean(np.abs(mcft_out_vec),3),2)
# mcft2_mean_vec = np.mean(np.mean(np.abs(mcft_out_vec)**2,3),2)
#
#
# # mcft_mean/= mcft_mean.max()
# # mcft_mean_vec/=mcft_mean_vec.max()
#
# print(mcft_mean.min(),mcft_mean.max())
# print(mcft2_mean.min(),mcft2_mean.max())
#
# print(mcft_mean_vec.min(),mcft_mean_vec.max())
# print(mcft2_mean_vec.min(),mcft2_mean_vec.max())


# plt.figure()
# plt.subplot(221)
# plt.pcolormesh(mcft_mean)
# plt.title('mcft_ft')
# plt.subplot(222)
# plt.pcolormesh(mcft2_mean)
# plt.title('mcft2_ft')
# plt.subplot(223)
# plt.pcolormesh(mcft_mean_vec)
# plt.title('mcft_sr')
# plt.subplot(224)
# plt.pcolormesh(mcft2_mean_vec)
# plt.title('mcft2_sr')
# plt.show()


################################## FFT speedup ##################################

# d0,d1,d2,d3 = np.shape(fbank_sr_domain)
# #nn = 10000
# aa = np.random.rand(d0,d1,d2,d3)
#
# nfft2 = helper.next_fast_len(d2)
# nfft3 = helper.next_fast_len(d3)
# print('nfft2: ',nfft2, ' nfft3: ',nfft3)
#
# bb1 = np.zeros((d0,d1,d2,d3), dtype='complex128')
# start = time.time()
# for i in range(d0):
#     for j in range(d1):
#         aa_temp = aa[i,j,:,:]
#         bb1_temp = np.fft.fft(aa_temp,axis=1)
#         bb1_temp = np.fft.fft(bb1_temp.T,axis=1)
#         bb1[i,j,:,:] = bb1_temp.T
#
# stop = time.time()
# print('computation time: ',stop-start)
#
# bb2 = np.zeros((d0,d1,d2,d3), dtype='complex128')
# start = time.time()
# for i in range(d0):
#     for j in range(d1):
#         aa_temp = aa[i, j, :, :]
#         bb2_temp = np.fft.fft(aa_temp, axis=1)
#         bb2_temp = np.fft.fft(bb2_temp.T, axis=1)
#         bb2[i, j, :, :] = bb2_temp.T
#
# stop = time.time()
# print('computation time: ',stop-start)
#
#
# err = np.sum(np.abs(bb1-bb2))
# print('err: ',err)