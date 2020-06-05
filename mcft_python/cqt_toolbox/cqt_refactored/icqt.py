import numpy as np
from filterbank import gen_inv_filterbank

def apply_inv_filterbank(cqt_in, inv_filter_bank, shift_ctr_freq, sig_len, phasemode):
    # type: (list, list, np.ndarray, int, str) -> np.ndarray
    '''
    Input parameters:
          cqt_in            : List (or ndarray) containing the transform coefficients (nonstationary Gabor coefficients)
          inv_filterbank    : List containing the frequency-doamin synthesis filters
          shift_ctr_freq    : Ndarray of center frequency shifts
          sig_len           : Length of the analyzed signal
          phasemode         : Can be set to (default is 'global')
                              - 'local'  : Zero-centered filtered used
                              - 'global' : Mapping function used (see reference)
    Output parameters:
          rec_signal        : Synthesized real-valued signal (Channels are stored
                              in the columns)

    This function synthesizes a real-valued signal by applying the provided nonstationary Gabor filterbank
    to the cqt coefficients positioned by the ndarray of shifts_ctr_freq.

    Note that, due to the structure of the coefficient array in the real valued setting, all entries g{n}
    with N > length(c) will be ignored and assumed to be fully supported on the negative frequencies.

    Let P(n)=sum_{l=1}^{n} shift(l), then the synthesis formula reads:

                     N-1
        fr_temp(l) = sum sum c{n}(m)g{n}[l-P(n)]*exp(-2*pi*i*(l-P(n))*m/M(n)),
                     n=0  m

    for l=0,...,Ls-1.  In practice, the synthesis formula is realized by fft and overlap-add. To synthesize
    the negative frequencies, fr_temp is truncated to length floor( Ls/2 )+1. Afterwards ifftreal implicitly
    computes the hermite symmetric extension and computes the inverse Fourier transform, i.e. fr = ifftreal(fr_temp).

    If a nonstationary Gabor frame was used to produce the coefficients and inv_filterbank is a corresponding
    dual frame, this function should perfectly reconstruct the originally analyzed signal to numerical precision.

    Multichannel output will save each channel in a column of rec_signal, such that the output is of the shape
    (timeseries, num_channels).

    References:
      C. Schorkhuber, A. Klapuri, N. Holighaus, and M. Dorfler. A Matlab
      Toolbox for Efficient Perfect Reconstruction log-f Time-Frequecy
      Transforms.

      P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
      Theory, implementation and applications of nonstationary Gabor Frames.
      J. Comput. Appl. Math., 236(6):1481-1496, 2011.

      G. A. Velasco, N. Holighaus, M. DAJrfler, and T. Grill. Constructing an
      invertible constant-Q transform with non-stationary Gabor frames.
      Proceedings of DAFX11, Paris, 2011.

    See also:  gen_inv_filterbank, icqt

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    '''
    # Input checking
    if len(cqt_in[0].shape) > 1:
        num_channels = cqt_in[0].shape[1]
    else:
        num_channels = 1
    num_filters = len(cqt_in)

    # Variable initialization for reconstructed signal
    ctr_freqs = np.cumsum(shift_ctr_freq)
    rec_sig_len = ctr_freqs[-1]
    ctr_freqs -= shift_ctr_freq[0]
    rec_signal = np.zeros((int(rec_sig_len), int(num_channels)), dtype=np.complex128)

    # Apply the inverse filterbank to generate the reconstructed signal
    for i in range(num_filters):
        filter_len = len(inv_filter_bank[i])
        filter_range = ((ctr_freqs[i] + np.arange(-1 * np.floor(filter_len / 2),
                                                  np.ceil(filter_len / 2))) % rec_sig_len).astype(np.int32)

        filt_sig_temp = np.fft.fft(cqt_in[i], axis=0) * len(cqt_in[i])

        if phasemode == 'global':
            fsNewBins = len(cqt_in[i])
            fkBins = ctr_freqs[i]
            displace = fkBins - np.floor(fkBins / fsNewBins) * fsNewBins
            filt_sig_temp = np.roll(filt_sig_temp, -1 * int(displace))

        # Apply the filter over the computed filter support (filter_range)
        first_half = np.arange(len(filt_sig_temp) - np.floor(filter_len / 2), len(filt_sig_temp))
        second_half = np.arange(np.ceil(filter_len / 2))
        idx1 = np.concatenate((first_half, second_half))
        filt_sig_temp = filt_sig_temp[np.mod(idx1, len(filt_sig_temp)).astype(np.int32)]
        idx2 = np.concatenate(
            (np.arange(filter_len - np.floor(filter_len / 2), filter_len), np.arange(np.ceil(filter_len / 2)))).astype(
            np.int32)
        rec_signal[filter_range, :] += np.expand_dims(filt_sig_temp * inv_filter_bank[i][idx2],axis=1)

    # Compute the " > nyq_bin " side of the spectrum using Hermitian symmetry
    nyq_bin = int(np.floor(sig_len / 2))
    rec_signal[nyq_bin + 1:] = np.conj(rec_signal[nyq_bin - (1 - (sig_len % 2)):0:-1])

    # Compute the inverse Fourier transform of the whole spectrum
    rec_signal = np.real(np.fft.ifft(rec_signal, axis=0))

    return rec_signal

def icqt(cqt_main, cqt_dc, cqt_nyq, filter_bank, fbank_params):
    # type: (np.ndarray,np.ndarray,np.ndarray,list,dict) -> (np.ndarray, list)
    '''
    Input parameters:
          cqt_main         : Numpy array containing cqt coefficients generated by applying a constant-Q or linear-Q
                             bandpass filterbank to the input time-domain signal.
          cqt_dc           : Numpy array containing the lowpass filtered part of the time signal.
          cqt_nyq          : Numpy array containing the highpass filtered part of the time signal.
          filter_bank      : List containing all filters (lowpass, highpass, bandpass over positive and negative frequencies)

          fbank_params              : Dictionary containing parameters required for the inverse transform, including
              .shift_ctr_freq     : Ndarray of center frequency shifts
              .bw_samp_int        : Ndarray containing filter bandwidths (in number of frequency samples)
              .sig_len            : Length of input signal
              .phasemode          : 'local'  --> zero-centered filtered used
                                      'global' --> frequency mapping function used
              .rasterize          : time-frequency-domain sampling scheme (full,piecewise,none)
              .bpass_ctr_freqs    : Numpy array containing the center frequencies of bandpass filters.
                                      It can also be used as the frequency vector for plotting purposes.

    Output parameters:
          rec_signal         : Ndarray containing the reconstructed time domain signal
          inv_filter_bank    : List containing the synthesis filterbank

    References:
      C. Schorkhuber, A. Klapuri, N. Holighaus, and M. Dorfler. A Matlab
      Toolbox for Efficient Perfect Reconstruction log-f Time-Frequecy
      Transforms.

      G. A. Velasco, N. Holighaus, M. Dorfler, and T. Grill. Constructing an
      invertible constant-Q transform with non-stationary Gabor frames.
      Proceedings of DAFX11, Paris, 2011.

      N. Holighaus, M. Dorfler, G. Velasco, and T. Grill. A framework for
      invertible, real-time constant-q transforms. Audio, Speech, and
      Language Processing, IEEE Transactions on, 21(4):775-785, April 2013.

    See also:  gen_inv_filterbank, apply_inv_filterbank

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    '''
    inv_filter_bank = gen_inv_filterbank(filter_bank, fbank_params['shift_ctr_freq'], fbank_params['bw_samp_int'])

    # We currently assume rasterize is always full
    cqt_in = [x for x in cqt_main]
    cqt_in.insert(0, cqt_dc)
    cqt_in.append(cqt_nyq)

    rec_signal = apply_inv_filterbank(cqt_in, inv_filter_bank, fbank_params['shift_ctr_freq'], fbank_params['sig_len'],
                                      fbank_params['phasemode'])

    return rec_signal, inv_filter_bank