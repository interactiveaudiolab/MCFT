from __future__ import print_function,division

import numpy as np


def gen_inv_filterbank(filter_bank,shift,bw_bins=None):
    # type: (list, numpy.ndarray, numpy.ndarray) -> list
    '''
    Input parameters:
          filter_bank       : List of filters
          shift             : Ndarray of time/frequency shifts
          bw_bins           : Number/Vector of frequency channels
    Output parameters:
          inv_filter_ban    : Canonical dual frame filters

    gen_inv_filterbank uses the forward filterbank and shift vector to compute
    the canonical dual frame windows/filters by inverting the diagonal of the 
    frame operator and applying the inverse to the forward filters. More explicitly, 

       gd{n} = g{n} / ( sum M(l) |g{l}|^2 ), 
                         l   

    If filterbank, shift, and bw_bins specify a painless frame, i.e. 

        SUPP(G{N})  <= M(n)~forall~n

    and

        A <= sum ( M(n) |g{n}|^2 ) <= B, for some 0 < A <= B < infty
              n   

    the computation will result in the canonical dual frame. If filterbank, 
    shift,and bw_bins specify a frame, but the first condition is violated,
    the result can be interpreted as a first approximation of the 
    corresponding canonical dual frame.
  
    Note, the time shifts corresponding to the dual window sequence is the
    same as the original shift sequence and as such already given. 

    If filterbank, shift, and bw_bins is a painless frame, the output can 
    be used for perfect reconstruction of a signal using the inverse 
    nonstationary Gabor transform in apply_inv_filterbank.

    References:
      P. Balazs, M. DAJrfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
      Theory, implementation and applications of nonstationary Gabor Frames.
      J. Comput. Appl. Math., 236(6):1481-1496, 2011.
  
    See also:  apply_inv_filterbank, icqt

    Translation from MATLAB by: Trent Cwiok (cwiok@u.northwestern.edu)
                                Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
    '''
    # Argument handling
    if bw_bins is None:
        bw_bins = np.zeros(len(shift))
        for i in range(len(shift)):
            bw_bins[i] = len(filter_bank[i])

    if bw_bins.size == 1:
        bw_bins = bw_bins[0]*np.ones(len(shift))

    # Initialize variables
    num_filters = len(shift)
    ctr_freqs = np.cumsum(shift)
    sig_len = ctr_freqs[num_filters-1]
    ctr_freqs -= shift[0]
    diag = np.zeros(int(sig_len))
    filter_range = []

    # Create the diagonal of the frame operator
    for i in range(num_filters):
        filter_len = len(filter_bank[i])
        filter_range.append(((ctr_freqs[i] + np.arange(-1*np.floor(filter_len/2),np.ceil(filter_len/2))) % sig_len).astype(np.int32))
        diag[filter_range[i]] += (np.fft.fftshift(filter_bank[i])**2)*bw_bins[i]

    # Compute the inverse filters and store them
    inv_filter_bank = filter_bank
    for i in range(num_filters):
        inv_filter_bank[i] = np.fft.ifftshift(np.fft.fftshift(inv_filter_bank[i])/diag[filter_range[i]])
    
    return inv_filter_bank