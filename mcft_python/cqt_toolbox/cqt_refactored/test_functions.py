import numpy as np
from filterbank import gen_filter, gen_filterbank, gen_inv_filterbank
from cqt import apply_filterbank, cqt
from icqt import apply_inv_filterbank, icqt

# generate a time signal
samp_rate = 8000
sig_len = 2 * samp_rate
tvec = np.arange(sig_len)/samp_rate
f1 = 27.5 * 2**(2/12)
signal = np.cos(2 * np.pi * f1 * tvec)
signal = np.expand_dims(signal,1)

# filterbank parameters
fmin = 27.5 * 2**(0/12)
fmax = 27.5 * 2**(12/12)
bins_per_octave = 12
phasemode = 'global'


# 1. gen_filter: generate the values of a single frequency-domain filter
filt = gen_filter('hann', num_samples=6)
print('filter values: ', filt)

# 2. gen_filterbank: generate a filterbank in the frequency domain
fbank, shift_ctr_freq, bw_samp_int = gen_filterbank(fmin, fmax, bins_per_octave, samp_rate, sig_len,
                                                          min_filt_len=4, bw_factor=1, fractional=False,
                                                          filt_type='hann', gamma=0)
# 3. apply_filterbank
cqt_out, sig_len = apply_filterbank(signal, fbank, shift_ctr_freq, phasemode, bw_samp_int=bw_samp_int)

# 4. gen_inv_filterbank
inv_fbank = gen_inv_filterbank(fbank,shift_ctr_freq,bw_samp_int=bw_samp_int)

# 5. apply_inv_filterbank
cqt_in = [np.squeeze(i) for i in cqt_out]
rec_signal_1 = apply_inv_filterbank(cqt_in, inv_fbank, shift_ctr_freq, sig_len, phasemode)

# 6. cqt
cqt_main, cqt_dc, cqt_nyq, fbank, fbank_params = cqt(signal, bins_per_octave, samp_rate, fmin, fmax, rasterize='full',
                                          phasemode='global', gamma=0, normalize='sine',filt_type='hann')

# 7. icqt
rec_signal_2, inv_filter_bank = icqt(cqt_main, cqt_dc, cqt_nyq, fbank, fbank_params)

# Compute reconstruction errors
def norm2(x):
    y = np.sqrt(np.sum(x**2))
    return y

rec_err_1 = 20 * np.log10(norm2(signal - rec_signal_1)/norm2(signal))
rec_err_2 = 20 * np.log10(norm2(signal - rec_signal_2)/norm2(signal))

print(rec_err_1, rec_err_2)