import numpy as np
from filterbank import gen_filterbank

def apply_filterbank(signal, filter_bank, shift_ctr_freq, phasemode, bw_samp_int=None):
    # type: (np.ndarray, list, np.ndarray, str, np.ndarray) -> (list, int)
    '''
    Input parameters:
        signal          : A real-valued signal -- multichannel
                          signals should be ndarray of shape
                          (len_of_signal, num_channels)
        filter_bank     : List containing the frequency-domain analysis filterbank (logarithmically spaced)
        shift_ctr_freq  : Ndarray of center frequency shifts
        phasemode       : 'local': zero-centered filtered used
                          'global': frequency mapping function used (see cqt)
        **bw_samp_int   : Number of time channels
                          If this is constant, the output is converted to a matrix

    Output parameters:
        cqt_out         : List (or ndarray) containing the transform coefficients
        sig_len         : Original signal length

    **optional args

    Given a (multichannel) signal, bank of filters, and frequency shifts, apply_filterbank applies
    the corresponding nonstationary Gabor filter to the real signal, using only the filters that are
    at least partially supported on the positive frequencies. Let P(n) = sum_{l=1}^{n} shift(l), then
    the output cqt = apply_filterbank(signal,filter_bank,shift,bw_bins) is a list with

              Ls-1
      c{n}(m)= sum fft(f)(l)*conj(g{n}(l-P(n)))*exp(2*pi*i*(l-P(n))*m/M(n))
               l=0

    where m runs from 0 to M(n)-1 and n from 1 to N, where g{N} is the final filter at least partially
    supported on the positive frequencies. All filters in the bank and shifts that are completely
    supported on the negative frequencies are ignored.

    See also:  gen_filterbank, gen_inv_filterbank, apply_inv_filterbank

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

    '''

    # Unpack the signal length and num of channels
    sig_len, num_channels = signal.shape

    # Generate some variables for computation later on
    num_filters = len(shift_ctr_freq)
    if bw_samp_int is None:
        bw_samp_int = np.zeros(num_filters)
        for i in range(num_filters):
            bw_samp_int[i] = len(filter_bank[i])

    if bw_samp_int.size == 1:
        bw_samp_int = bw_samp_int[0] * np.ones(num_filters)

    # Convert the time-domain signal to frequency-domain
    signal = np.fft.fft(signal, axis=0)

    # Compute the center frequencies of filters using the shifts between center frequencies
    ctr_freqs = np.cumsum(shift_ctr_freq) - shift_ctr_freq[0]

    # A small amount of zero-padding might be needed
    zpad_len = np.sum(shift_ctr_freq) - sig_len
    padding = np.zeros((int(zpad_len), int(num_channels)))
    signal = np.vstack((signal, padding))

    # Compute the filter lengths
    filter_lens = np.zeros(len(filter_bank))
    for i in range(len(filter_bank)):
        filter_lens[i] = len(filter_bank[i])

    # Number of filters determined by the last filter at least partially supported on frequencies
    # below the Nyquist rate (find the highest filter whose lower end is smaller that the Nyquist rate)
    num_filters = [ctr_freqs[i] - np.floor(filter_lens[i] / 2) <= (sig_len + zpad_len) / 2 for i in
                   range(len(ctr_freqs))]
    num_filters = np.nonzero(num_filters)[0][-1]

    # Applying the filters at the correct center frequencies
    cqt_out = []
    for i in range(num_filters + 1):
        # Create the lists of indices at which points the filters will be applied to the signal
        idx = np.concatenate(
            (np.arange(np.ceil(filter_lens[i] / 2), filter_lens[i]), np.arange(np.ceil(filter_lens[i] / 2))))
        filter_range = ((ctr_freqs[i] + np.arange(-1 * np.floor(filter_lens[i] / 2),
                                                  np.ceil(filter_lens[i] / 2))) % sig_len + zpad_len)
        idx, filter_range = (idx.astype(np.int32), filter_range.astype(np.int32))

        # Create the cqt coefficients
        if bw_samp_int[i] < filter_lens[i]:
            # This case involves aliasing (non-painless case)
            num_columns = np.ceil(filter_lens[i] / bw_samp_int[i])
            filt_out_temp = np.zeros((int(num_columns) * int(bw_samp_int[i]), num_channels), dtype=np.complex128)

            idx_slice_one = np.arange((filt_out_temp.shape[0] - np.floor(filter_lens[i] / 2)), filt_out_temp.shape[0], dtype=np.int32)
            idx_slice_two = np.arange(np.ceil(filter_lens[i] / 2), dtype=np.int32)

            filt_out_temp[np.concatenate((idx_slice_one, idx_slice_two)), :] = signal[filter_range, :] * filter_bank[i][idx]

            filt_out_temp = np.reshape(filt_out_temp, (bw_samp_int[i], num_columns, num_channels))
            cqt_out.append(np.squeeze(np.fft.ifft(np.sum(filt_out_temp, axis=1))))

        else:
            # No aliasing here
            filt_out_temp = np.zeros((int(bw_samp_int[i]), num_channels), dtype=np.complex128)
            idx_slice_one = np.arange((filt_out_temp.shape[0] - np.floor(filter_lens[i] / 2)), filt_out_temp.shape[0], dtype=np.int32)
            idx_slice_two = np.arange(np.ceil(filter_lens[i] / 2), dtype=np.int32)
            filt_out_temp[np.concatenate((idx_slice_one, idx_slice_two)), :] = signal[filter_range, :] * \
                                                                               np.reshape(filter_bank[i][idx],(len(filter_bank[i]),1))

            if phasemode == 'global':
                fsNewBins = bw_samp_int[i]
                fkBins = ctr_freqs[i]
                displace = fkBins - np.floor(fkBins / fsNewBins) * fsNewBins
                filt_out_temp = np.roll(filt_out_temp, int(displace))

            cqt_out.append(np.fft.ifft(filt_out_temp, axis=0))

    # If the outputs of all filters are of the same length, reshape the list into an ndarray
    if np.max(bw_samp_int) == np.min(bw_samp_int):
        cqt_out = np.asarray(cqt_out)
        cqt_out = np.reshape(cqt_out, (int(bw_samp_int[0]), int(num_filters + 1), int(num_channels)))

    return cqt_out, sig_len


def cqt(signal, bins_per_octave, samp_rate, fmin, fmax, rasterize='full', phasemode='global', gamma=0, normalize='sine',
        filt_type='hann'):
    # type: (np.ndarray, int, int, float, float, str, str, int, str, str) -> (np.ndarray,np.ndarray,np.ndarray,list,dict)
    '''
    Input parameters:
          signal            : A real-valued signal -- multichannel signals should be ndarray of shape
                              (len_of_signal, num_channels)
          bins_per_octave   : Number of bandpass filters per octave
          samp_rate         : Sampling frequency
          fmin              : Lowest frequency to be analyzed
          fmax              : Highest frequency to be analyzed
          **rasterize       : Can be 'none', 'full', or 'piecewise' -- affects the hop size (subsampling rate)
                              Note: only 'full' is implemented in the current version
                              'full' --> The hop sizes for all freqency channels are set to the smallest
                                         hop size in the representation. Transform coefficients will be
                                         presented in matrix format.
          **phasemode       : Can be 'local' or 'global' -- default: 'global'
                              'local' --> Zero-centered filtered used
                              'global' --> Mapping function used (see reference)
          **gamma           : Constant factor for offsetting filter bandwidths, >=0
                              The bandwidth of each filter is given by Bk = 1/Q * fk + gamma, where
                              fk is the filters center frequency, Q is fully determined by the number
                              of bins per octave and gamma is a bandwidth offset. If gamma = 0 the obtained
                              filterbank is constant-Q. Setting gamma > 0 time resolution towards lower
                              frequencies can be improved compared to the constant-Q case
                              (e.g. ERB proportional bandwidths). See reference for more information.
          **normalize       : Can be 'sine', 'impulse' or 'none' -- used to normalize the coefficients
                              'sine'    --> Filters are scaled such that a sinusoid with amplitude A in time domain
                                            will exhibit the same amplitude in the time-frequency representation.
                              'impulse' --> Filters are scaled such that an impulse in time domain will exhibit
                                            a flat response in the time-frequency representation (in the frame
                                            that centers the impulse)
          **filt_type      : String specifying the name of the weighting function to be used as a filter.
                             It can be a real-valued window, e.g., 'hann'. See gen_filter docstring for a
                             list of available filter types.

    Output parameters:
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

    **optional args

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

    See also:  gen_filter, gen_filterbank, apply_filterbank

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    '''

    # Generate the frequency-domain analysis filterbank
    filter_bank, shift_ctr_freq, bw_samp_int = gen_filterbank(fmin, fmax, bins_per_octave, samp_rate, len(signal),
                                                 filt_type=filt_type, gamma=gamma)

    # Compute the number and center frequencies of bandpass filters
    num_bpass_filters = int(len(bw_samp_int) / 2 - 1)
    bpass_ctr_freqs = samp_rate * np.cumsum(shift_ctr_freq[1:]) / len(signal)
    bpass_ctr_freqs = bpass_ctr_freqs[:num_bpass_filters]

    # Assume all bandpass filters have the same number of samples (rasterize = 'full')
    bw_samp_int[1:num_bpass_filters + 1] = bw_samp_int[num_bpass_filters]
    bw_samp_int[num_bpass_filters + 2:] = bw_samp_int[num_bpass_filters:0:-1]

    # Create a normalization vector
    normalize = normalize.lower()
    if normalize in ['sine', 'sin']:
        norm_vec = 2 * bw_samp_int[:num_bpass_filters + 2] / len(signal)
    elif normalize in ['impulse', 'imp']:
        filter_lens = np.array([len(i) for i in filter_bank])
        norm_vec = 2 * bw_samp_int[:num_bpass_filters + 2] / filter_lens[:num_bpass_filters + 2]
    elif normalize in ['none', 'no']:
        norm_vec = np.ones(num_bpass_filters + 2)
    else:
        print("Unknown normalization method")

    norm_vec = np.concatenate((norm_vec, norm_vec[len(norm_vec) - 2:0:-1]))

    # Apply normalization to the filterbank
    for i in range(len(norm_vec)):
        filter_bank[i] *= norm_vec[i]

    if len(signal.shape) < 2:
        signal.shape = (len(signal), 1)

    # Apply the normalized filterbank to the signal to compute the cqt
    cqt, sig_len = apply_filterbank(signal, filter_bank, shift_ctr_freq, phasemode, bw_samp_int)

    # Assume rasterize is full always
    # Store the lowpass, highpass, and bandpass filtered sections separately (cqt matrix is the bandpass section)
    cqt_dc = np.squeeze(cqt[0])
    cqt_nyq = np.squeeze(cqt[num_bpass_filters + 1])
    cqt_main = np.squeeze(np.asarray(cqt[1:num_bpass_filters + 1]))

    fbank_params = {'shift_ctr_freq': shift_ctr_freq, 'bw_samp_int': bw_samp_int, 'sig_len': sig_len,
               'phasemode': phasemode, 'rasterize': rasterize, 'bpass_ctr_freqs': bpass_ctr_freqs}

    return cqt_main, cqt_dc, cqt_nyq, filter_bank, fbank_params