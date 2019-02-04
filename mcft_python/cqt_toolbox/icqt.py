from __future__ import print_function,division

import numpy as np

from gen_inv_filterbank import gen_inv_filterbank
from apply_inv_filterbank import  apply_inv_filterbank


def icqt(Xcq):
    # type: (dict) -> (ndarray, list)
    '''
    Input parameters:
          Xcq                : dict obtained by cqt(...)

    Output parameters: 
          rec_signal         : reconstructed time domain signal
          inv_filter_bank    : synthesis filterbank 

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
    Xcq['inv_filter_bank'] = gen_inv_filterbank(Xcq['filter_bank'],Xcq['shift'],Xcq['bw_bins'])

    # We currently assume rasterize is always full
    cqt = [x for x in Xcq['cqt']]
    cqt.insert(0,Xcq['cqt_DC'])
    cqt.append(Xcq['cqt_Nyq'])
    
    rec_signal = apply_inv_filterbank(cqt,Xcq['inv_filter_bank'],Xcq['shift'],Xcq['sig_len'],Xcq['phasemode'])
    
    inv_filter_bank = Xcq['inv_filter_bank']
    
    return rec_signal, inv_filter_bank