import numpy as np
from scipy.fftpack import fft2,ifft2
from cqt_toolbox.icqt import icqt

def inv_mcft(mcft_in, cqt_params, fbank_sr_domain):
    """
    This function reconstructs a time-doman signal given its
    Multi-resolution Common Fate Transform (MCFT).
    The intermediary time-frequency domain representation is the
    Constant-Q Transform (CQT) of the audio signal, which is computed
    using the invertible and optimized CQT implementation proposed by
    Schorkhuber et al.:

    Toolbox webpage:
    http://www.cs.tut.fi/sgn/arg/CQT/
    Reference:
    Schörkhuber et al. "A Matlab toolbox for efficient perfect reconstruction
    time-frequency transforms with log-frequency resolution."
    Audio Engineering Society Conference: 53rd International Conference:
    Semantic Audio. Audio Engineering Society, 2014.

    The python translation of this CQT toolbox can be found here:
    https://github.com/interactiveaudiolab/MCFT/tree/master/mcft_python/cqt_toolbox

    Inputs:
    mcft_in: 4d numpy array containing the MCFT coefficients
    cqt_params: dictionary containing cqt parameters, including:
                num_freq_bin: total number of frequency bins
                num_time_frame: total number of time frames

    fbank_sr_domain: 4d numpy array containing the scale-rate-domain filterbank

    Output:
    signal_rec: numpy array containing the reconstructed time-doamin signal

    Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    """

    ### reconstruct the signal cqt

    print('Reconstructing the CQT...')
    sig_cqt_rec = mcft_to_cqt(mcft_in,fbank_sr_domain)

    ### reconstruct the time-doamin signal

    print('Reconstructing the time-domain signal...')

    num_freq_bin = cqt_params['num_freq_bin']
    num_time_frame = cqt_params['num_time_frame']

    del cqt_params['num_freq_bin']
    del cqt_params['num_time_frame']

    cqt_dict_full = cqt_params
    cqt_dict_full['cqt'] = sig_cqt_rec[0:num_freq_bin, 0:num_time_frame]

    signal_rec, _  = icqt(cqt_dict_full)

    return signal_rec



def mcft_to_cqt(mcft_in,fbank_sr_domain):
    """
    This function reconstructs the time-frequency representation (CQT)
    of an audio signal through inverse filtering given the 4-dimensional
    MCFT representation and the scale-rate-domain filterbank.

    Inputs:
    mcft_in: 4d numpy array containing the MCFT coefficients
    fbank_sr_domain: 4d numpy array containing a bank of scale-rate-domain filters

    Ouput:
    sig_cqt_rec: 2d numpy array containing the reconstructed time-frequency
                 representation (complex in general)


    Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    """

    # dimensions
    num_scale_ctrs, num_rate_ctrs, nfft_scale, nfft_rate = np.shape(mcft_in)

    # MCFT to time-frequency representation

    fbank_sr_sum = np.zeros((nfft_scale,nfft_rate),dtype='complex128')
    cqt_2dft_sum = np.zeros((nfft_scale,nfft_rate),dtype='complex128')
    for i in range(num_scale_ctrs):
        for j in range(num_rate_ctrs):

            mcft_temp = mcft_in[i,j,:,:]
            cqt_2dft_temp = fft2(mcft_temp,[nfft_scale,nfft_rate])

            filt_sr_temp = fbank_sr_domain[i,j,:,:]

            cqt_inv_filt = cqt_2dft_temp * np.conj(filt_sr_temp)

            fbank_sr_sum += filt_sr_temp * np.conj(filt_sr_temp)
            cqt_2dft_sum += cqt_inv_filt


    # normalize cqt_2dft_sum by fbank_sr_sum
    cqt_2dft_ratio = cqt_2dft_sum / (fbank_sr_sum+1e-16)

    # compute the reconstructed cqt
    sig_cqt_rec = ifft2(cqt_2dft_ratio)

    return sig_cqt_rec


