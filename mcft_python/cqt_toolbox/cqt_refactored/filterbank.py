import numpy as np

def gen_filter(filt_type, sample_positions=None, num_samples=None):
    # type: (str, np.ndarray, float) -> np.ndarray
    '''
    Input parameters:
          filt_type          : String specifying the name of the weighting function to be used as a filter.
                               It can be a real-valued window, e.g., 'hann'.
          **sample_positions   : Ndarray containing the positions of filter samples in the frequency domain.
          **num_samples        : Filter length (in number of samples).

    Output parameters:
          filter_              : Ndarray containing the frequency-domain filter.

    **optional args

    This function is used to generate an individual filter in the shape of the specified window.
    If sample_positions is NOT None, the filter is calculated at the specified sample positions.
    If sample_positions is None, then num_samples must be provided. In this case, the filter is calculated
    over the interval [-0.5,0.5) at num_samples equally spaced sampling positions.
    Before returning the filter it is masked to force everything outside of -.5 and .5 to zero.

    The following filter types are available:

    'hann'       von Hann window. Forms a PU. The Hann window has a
                 mainlobe with of 8/N, a PSL of -31.5 dB and decay rate
                 of 18 dB/Octave.

    'cos'        Cosine window. This is the square root of the Hanning
                 window. The cosine window has a mainlobe width of 6/N,
                 a  PSL of -22.3 dB and decay rate of 12 dB/Octave.

    'rec'        Rectangular window. The rectangular window has a
                 mainlobe width of 4/N, a  PSL of -13.3 dB and decay
                 rate of 6 dB/Octave. Forms a PU. Alias: 'square'

    'tri'        Triangular window.

    'hamming'    Hamming window. Forms a PU that sums to 1.08 instead
                 of 1.0 as usual. The Hamming window has a
                 mainlobe width of 8/N, a  PSL of -42.7 dB and decay
                 rate of 6 dB/Octave.

    'blackman'   Blackman window. The Blackman window has a
                 mainlobe width of 12/N, a PSL of -58.1 dB and decay
                 rate of 18 dB/Octave.

    'blackharr'  Blackman-Harris window. The Blackman-Harris window has
                 a mainlobe width of 16/N, a PSL of -92.04 dB and decay
                 rate of 6 dB/Octave.

    'modblackharr'  Modified Blackman-Harris window. This slightly
                    modified version of the Blackman-Harris window has
                    a mainlobe width of 16/N, a PSL of -90.24 dB and decay
                    rate of 18 dB/Octave.

    'nuttall'      Nuttall window. The Nuttall window has a mainlobe
                   width of 16/N, a PSL of -93.32 dB and decay rate of
                   18 dB/Octave.

    'nuttall10'    2-term Nuttall window with 1 continuous derivative.
                   Alias: 'hann'.

    'nuttall01'    2-term Nuttall window with 0 continuous derivatives.
                   Alias: 'hamming'.

    'nuttall20'    3-term Nuttall window with 3 continuous derivatives.
                   The window has a mainlobe width of 12/N, a PSL of
                   -46.74 dB and decay rate of 30 dB/Octave.

    'nuttall11'    3-term Nuttall window with 1 continuous derivative.
                   The window has a mainlobe width of 12/N, a PSL of
                   -64.19 dB and decay rate of 18 dB/Octave.

    'nuttall02'    3-term Nuttall window with 0 continuous derivatives.
                   The window has a mainlobe width of 12/N, a PSL of
                   -71.48 dB and decay rate of 6 dB/Octave.

    'nuttall30'    4-term Nuttall window with 5 continuous derivatives.
                   The window has a mainlobe width of 16/N, a PSL of
                   -60.95 dB and decay rate of 42 dB/Octave.

    'nuttall21'    4-term Nuttall window with 3 continuous derivatives.
                   The window has a mainlobe width of 16/N, a PSL of
                   -82.60 dB and decay rate of 30 dB/Octave.

    'nuttall12'    4-term Nuttall window with 1 continuous derivatives.
                   Alias: 'nuttall'.

    'nuttall03'    4-term Nuttall window with 0 continuous derivatives.
                   The window has a mainlobe width of 16/N, a PSL of
                   -98.17 dB and decay rate of 6 dB/Octave.

    'gauss'        Truncated, stretched Gaussian: exp(-18*x^2) restricted
                   to the interval ]-.5,.5[.

    'wp2inp'       Warped Wavelet uncertainty equalizer (see WP 2 of the
                   EU funded project UnlocX). This function is included
                   as a test function for the Wavelet transform
                   implementation and serves no other purpose in this
                   toolbox.

    References:
      Wikipedia. Window function - wikipedia article.
      http://en.wikipedia.org/wiki/Window_function.

      A. Nuttall. Some windows with very good sidelobe behavior. IEEE Trans.
      Acoust. Speech Signal Process., 29(1):84-91, 1981.

      F. Harris. On the use of windows for harmonic analysis with the
      discrete Fourier transform. Proceedings of the IEEE, 66(1):51 - 83,
      January 1978.

    See also:  gen_filterbank, gen_inv_filterbank

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    '''

    # Input argument checking -- either sample_positions or num_samples must be provided
    if sample_positions is not None:
        pass

    elif num_samples is not None:
        step = 1 / num_samples
        # Two cases if the num of samples is even or odd
        if num_samples % 2 == 0:
            # For even N the sampling interval is [0,.5-1/N] + [-.5,0) --> in this exact order
            first_half = np.linspace(0, .5 - step, int(num_samples * .5))
            second_half = np.linspace(-.5, -step, int(num_samples * .5))
            sample_positions = np.concatenate((first_half, second_half))
        else:
            # For odd N the sampling interval is [0,.5-1/(2N)] + [-.5+1/(2N),0)
            first_half = np.linspace(0, .5 - .5 * step, int(num_samples * .5) + 1)
            second_half = np.linspace(-.5 + .5 * step, -step, int(num_samples * .5))
            sample_positions = np.concatenate((first_half, second_half))
    else:
        print("Error: either sample_positions or num_samples must be provided.")
        return None

    # Switch case for possible window names, forcing everything to lowercase
    window_name = filt_type.lower()
    if window_name in ['hann', 'nuttall10']:
        filter_ = .5 + .5 * np.cos(2 * np.pi * sample_positions)

    elif window_name in ['cosine', 'cos', 'sqrthann']:
        filter_ = np.cos(np.pi * sample_positions)

    elif window_name in ['hamming', 'nuttall01']:
        filter_ = .54 + .46 * np.cos(2 * np.pi * sample_positions)

    elif window_name in ['square', 'rec', 'boxcar']:
        filter_ = np.asarray([int(abs(i) < .5) for i in sample_positions], dtype=np.float64)

    elif window_name in ['tri', 'triangular', 'bartlett']:
        filter_ = 1 - 2 * abs(sample_positions)

    elif window_name in ['blackman']:
        filter_ = .42 + .5 * np.cos(2 * np.pi * sample_positions) + .08 * np.cos(4 * np.pi * sample_positions)

    elif window_name in ['blackharr']:
        filter_ = .35875 + .48829 * np.cos(2 * np.pi * sample_positions) + .14128 * np.cos(
            4 * np.pi * sample_positions) + \
                  .01168 * np.cos(6 * np.pi * sample_positions)

    elif window_name in ['modblackharr']:
        filter_ = .35872 + .48832 * np.cos(2 * np.pi * sample_positions) + .14128 * np.cos(
            4 * np.pi * sample_positions) + \
                  .01168 * np.cos(6 * np.pi * sample_positions)

    elif window_name in ['nuttall', 'nuttall12']:
        filter_ = .355768 + .487396 * np.cos(2 * np.pi * sample_positions) + .144232 * np.cos(
            4 * np.pi * sample_positions) + \
                  .012604 * np.cos(6 * np.pi * sample_positions)

    elif window_name in ['nuttall20']:
        filter_ = 3 / 8 + 4 / 8 * np.cos(2 * np.pi * sample_positions) + 1 / 8 * np.cos(4 * np.pi * sample_positions)

    elif window_name in ['nuttall11']:
        filter_ = .40897 + .5 * np.cos(2 * np.pi * sample_positions) + .09103 * np.cos(4 * np.pi * sample_positions)

    elif window_name in ['nuttall02']:
        filter_ = .4243801 + .4973406 * np.cos(2 * np.pi * sample_positions) + .0782793 * np.cos(
            4 * np.pi * sample_positions)

    elif window_name in ['nuttall30']:
        filter_ = 10 / 32 + 15 / 32 * np.cos(2 * np.pi * sample_positions) + 6 / 32 * np.cos(
            4 * np.pi * sample_positions) + \
                  1 / 32 * np.cos(6 * np.pi * sample_positions)

    elif window_name in ['nuttall21']:
        filter_ = .338946 + .481973 * np.cos(2 * np.pi * sample_positions) + .161054 * np.cos(
            4 * np.pi * sample_positions) + \
                  .018027 * np.cos(6 * np.pi * sample_positions)

    elif window_name in ['nuttall03']:
        filter_ = .3635819 + .4891775 * np.cos(2 * np.pi * sample_positions) + .1365995 * np.cos(
            4 * np.pi * sample_positions) + \
                  .0106411 * np.cos(6 * np.pi * sample_positions)

    elif window_name in ['gauss', 'truncgauss']:
        filter_ = np.exp(-18 * sample_positions ** 2)

    elif window_name in ['wp2inp']:
        filter_ = np.exp(np.exp(-2 * sample_positions) * 25. * (1 + 2 * sample_positions))
        filter_ = filter_ / np.max(filter_)
    else:
        print("Error: unknown filter type ", filt_type, " provided.")
        return None

    # List comprehension makes a mask to force values outside -.5 and .5 to zero
    mask = [int(abs(i) < .5) for i in sample_positions]
    filter_ *= mask

    return filter_



def gen_filterbank(fmin, fmax, bins_per_octave, samp_rate, sig_len, min_filt_len=4, bw_factor=1, fractional=False,
                   filt_type='hann', gamma=0):
    # type: (float, float, int, int, int, int, int, bool, str, int) -> (list, np.ndarray, np.ndarray)
    '''
    Input parameters:
        fmin                 : Center frequency of the lowest bandpass filter (in Hz)
        fmax                 : Center frequency of the highest bandpass filter (in Hz)
        bins_per_octave      : Number of bandpass filters per octave
        samp_rate            : Sampling rate (in Hz)
        sig_len              : Length of the time-doamin signal (in samples)
        **min_filt_len       : Minimum filter length allowed (in number of samples)
        **bw_factor          : Values in bw_bins are rounded to factors of this
        **fractional         : True if shifts (of filter centers) can have fractional values
        **filt_type          : String specifying the name of the weighting function to be used as a filter.
                               It can be a real-valued window, e.g., 'hann'.
        **gamma              : Constant factor for offsetting filter bandwidths, >=0

    Output parameters:
        filter_bank         : List of ndarrays, each containing a different filter
        shift_ctr_freq      : Ndarray of shifts between the center frequencies
        bw_samp_int         : Ndarray containing filter bandwidths (in number of frequency samples)

    **optional args

    Create a nonstationary Gabor filterbank with a constant or varying Q-factor and relevant
    frequency range from fmin to fmax. To allow for perfect reconstruction, the frequencies
    outside that range will be captured by 2 additional filters placed on the zero and Nyquist
    frequencies, respectively.

    The Q-factor (quality factor) is the ratio of center frequency to bandwidth -->  cent_freq/bandwidth.
    The Q-factor is determined only by the bins_per_octave and the value of gamma. A gamma value of 0 implies
    a constant-Q filter. For gamma greater than 0, variable-Q filters are returned, which allow for better
    time-resolution in lower frequencies.

    For more details on the construction of the constant-Q nonstationary Gabor filterbank, please see the references.

    References:
        C. Schorkhuber, A. Klapuri, N. Holighaus, and M. Dorfler. A Matlab
        Toolbox for Efficient Perfect Reconstruction log-f Time-Frequecy
        Transforms.

        G. A. Velasco, N. Holighaus, M. DAJrfler, and T. Grill. Constructing an
        invertible constant-Q transform with non-stationary Gabor frames.
        Proceedings of DAFX11, Paris, 2011.

        N. Holighaus, M. DAJrfler, G. Velasco, and T. Grill. A framework for
        invertible, real-time constant-q transforms. Audio, Speech, and
        Language Processing, IEEE Transactions on, 21(4):775-785, April 2013.

    See also: gen_filter, gen_inv_filterbank

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    '''

    # Make sure fmax does not exceed the Nyquist rate
    nyquist = samp_rate / 2
    if fmax > nyquist:
        fmax = nyquist

    # Calculate the resolution of the frequency grid
    fftres = samp_rate / sig_len

    # Calculate the total number of bandpass filters a.k.a frequency bins
    bins_total = np.floor(bins_per_octave * np.log2(fmax / fmin))

    # Calculate the center frequencies of bandpass filters
    ctr_freqs_hz = fmin * 2 ** (np.asarray(np.arange(bins_total + 1)) / bins_per_octave)

    # Calculate the bandwidths of bandpass filters
    Q = 2 ** (1 / bins_per_octave) - 2 ** (-1 / bins_per_octave)
    filt_bw_hz = Q * ctr_freqs_hz + gamma

    # Make sure the support of highest filter will not exceed the Nyquist rate
    nonzeroIndices = np.nonzero([int(ctr_freqs_hz[i] + filt_bw_hz[i] / 2 > nyquist) for i in range(len(ctr_freqs_hz))])
    if nonzeroIndices[0].size > 0:
        ctr_freqs_hz = ctr_freqs_hz[:nonzeroIndices[0][0]]
        filt_bw_hz = filt_bw_hz[:nonzeroIndices[0][0]]

    # Make sure the support of the lowest filter will not go below zero
    nonzeroIndices = np.nonzero([int(ctr_freqs_hz[i] - filt_bw_hz[i] / 2 < 0) for i in range(len(ctr_freqs_hz))])
    if nonzeroIndices[0].size > 0:
        ctr_freqs_hz = ctr_freqs_hz[nonzeroIndices[0][-1]:]
        filt_bw_hz = filt_bw_hz[nonzeroIndices[0][-1]:]
        print("Warning: fmin set to ", fftres * np.floor(ctr_freqs_hz[0]/fftres), " Hz!")

    # Get the number of bandpass filters
    num_ctr_freqs = len(ctr_freqs_hz)

    # Include the lowpass filter (centered at 0) and the highpass filter (centered at Nyquist) --> [0,pi]
    # Then include mirrored bandpass filters --> (pi,2pi)
    ctr_freqs_hz = np.concatenate(([0], ctr_freqs_hz, [nyquist], samp_rate - np.flip(ctr_freqs_hz, 0)))
    bw_hz = np.concatenate(
        ([2 * fmin], filt_bw_hz, [ctr_freqs_hz[num_ctr_freqs + 2] - ctr_freqs_hz[num_ctr_freqs]], np.flip(filt_bw_hz, 0)))

    # Convert bandwidths and center frequencies from Hz to number of frequency samples
    bw_samp = bw_hz / fftres
    ctr_freqs_samp = ctr_freqs_hz / fftres

    # Find center positions of filters on the DFT grid (frequencies of analysis)
    ctr_freqs_int = np.zeros(ctr_freqs_samp.size)
    ctr_freqs_int[:num_ctr_freqs + 2] = np.floor(ctr_freqs_samp[:num_ctr_freqs + 2])
    ctr_freqs_int[num_ctr_freqs + 2:] = np.ceil(ctr_freqs_samp[num_ctr_freqs + 2:])

    # Find the distance between center frequencies of the filters
    shift_ctr_freq = np.concatenate(([-1 * ctr_freqs_int[-1] % sig_len], np.diff(ctr_freqs_int)))

    # Calculate the fractional portion of center frequency values if the 'fractional' option is used
    if fractional:
        frac_shift = ctr_freqs_samp - ctr_freqs_int
        bw_samp_int = np.ceil(bw_samp + 1)
    else:
        bw_samp = np.round(bw_samp)
        bw_samp_int = bw_samp

    # Set minimum bandwidth
    for i in range(2 * (num_ctr_freqs + 1)):
        if bw_samp[i] < min_filt_len:
            bw_samp[i] = min_filt_len
            bw_samp_int[i] = bw_samp[i]

    # Generate the actual filterbank using the gen_filter function
    if fractional:
        filter_bank = []
        for i in range(len(bw_samp_int)):
            samples = np.concatenate(
                (range(int(np.ceil(bw_samp_int[i] / 2) + 1)),
                 range(int(-1 * np.floor(bw_samp_int[i] / 2)), 0))).astype(float)
            samples -= frac_shift[i]
            samples /= bw_samp[i]
            filt = gen_filter(filt_type, sample_positions=samples)
            filt /= np.sqrt(bw_samp[i])
            filter_bank.append(filt)
    else:
        filter_bank = []
        for i in range(len(bw_samp)):
            filter_bank.append(gen_filter(filt_type, num_samples=bw_samp[i]))

    # Round to multiples of the bw factor input arg
    bw_samp_int = bw_factor * np.ceil(bw_samp_int / bw_factor)

    # Generate samples of the lowpass and highpass filters, both defined by a Tukey window
    for i in [0, num_ctr_freqs + 1]:
        if bw_samp_int[i] > bw_samp_int[i + 1]:
            filter_bank[i] = np.ones(int(bw_samp_int[i]))
            start = int(np.floor(bw_samp_int[i] / 2) - np.floor(bw_samp_int[i + 1] / 2))
            end = int(np.floor(bw_samp_int[i] / 2) + np.ceil(bw_samp_int[i + 1] / 2))
            filter_bank[i][start:end] = gen_filter('hann', num_samples=bw_samp_int[i + 1])
            filter_bank[i] /= np.sqrt(bw_samp_int[i])

    return filter_bank, shift_ctr_freq, bw_samp_int


def gen_inv_filterbank(filter_bank, shift_ctr_freq, bw_samp_int=None):
    # type: (list, np.ndarray, np.ndarray) -> list
    '''
    Input parameters:
          filter_bank       : List of ndarrays, each containing a different filter
          shift_ctr_freq    : Ndarray of shifts between the center frequencies
          bw_samp_int       : Ndarray containing filter bandwidths (in number of frequency samples)

    Output parameters:
          inv_filter_bank   : Dual window functions used for inverse filtering

    gen_inv_filterbank uses the analysis filterbank and the vector of center frequency shifts to compute
    the canonical dual frame windows/filters by inverting the diagonal of the frame operator and applying
    the inverse to the analysis filters. More explicitly,

       gd{n} = g{n} / ( sum M(l) |g{l}|^2 ), --> gd: inv_filt and g:filt
                         l

    If filterbank, shift_ctr_freq, and bw_samp_int specify a painless frame, i.e.

        SUPP(G{N})  <= M(n)~forall~n

    and

        A <= sum ( M(n) |g{n}|^2 ) <= B, for some 0 < A <= B < infty
              n

    the computation will result in the canonical dual frame. If filterbank, shift_ctr_freq, and bw_samp_int
    specify a frame, but the first condition is violated, the result can be interpreted as a first approximation
    of the corresponding canonical dual frame.

    Note, the time shifts corresponding to the dual window sequence is the same as the original shift sequence
    and as such already given.

    If filterbank, shift_ctr_freq, and bw_samp_int is a painless frame, the output can be used for perfect
    reconstruction of a signal using the inverse nonstationary Gabor transform in apply_inv_filterbank.

    References:
      P. Balazs, M. DAJrfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
      Theory, implementation and applications of nonstationary Gabor Frames.
      J. Comput. Appl. Math., 236(6):1481-1496, 2011.

    See also:  apply_inv_filterbank, icqt

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    '''

    # Check the input arguments
    if bw_samp_int is None:
        bw_samp_int = np.zeros(len(shift_ctr_freq))
        for i in range(len(shift_ctr_freq)):
            bw_samp_int[i] = len(filter_bank[i])

    if bw_samp_int.size == 1:
        bw_samp_int = bw_samp_int[0] * np.ones(len(shift_ctr_freq))

    # Initialize variables
    num_filters = len(shift_ctr_freq)
    ctr_freqs = np.cumsum(shift_ctr_freq)
    sig_len = ctr_freqs[num_filters - 1]
    ctr_freqs -= shift_ctr_freq[0]
    diag = np.zeros(int(sig_len))
    filter_range = []

    # Create the diagonal of the frame operator
    for i in range(num_filters):
        filter_len = len(filter_bank[i])
        filter_range.append(
            ((ctr_freqs[i] + np.arange(-1 * np.floor(filter_len / 2), np.ceil(filter_len / 2))) % sig_len).astype(
                np.int32))
        diag[filter_range[i]] += (np.fft.fftshift(filter_bank[i]) ** 2) * bw_samp_int[i]

    # Compute the inverse filters and store them
    inv_filter_bank = filter_bank
    for i in range(num_filters):
        inv_filter_bank[i] = np.fft.ifftshift(np.fft.fftshift(inv_filter_bank[i]) / diag[filter_range[i]])

    return inv_filter_bank