from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt


def gen_filter(window_name, sample_positions=None, sample_num=None):
    # type: (str, numpy.ndarray, float) -> numpy.ndarray
    '''
    Input parameters: 
          window_name          : String containing the window name
          sample_positions     : Ndarray of sampling positions, optional
          sample_num           : Number of samples in the window, optional
    Output parameters:
          g                    : Ndarray filter in specified window shape
    
    This function is used to generate an individual filter in the shape of the
    specified window. The filter can either be sampled over the specified 
    vector of points if sample_positions is NOT None. Before returning the
    filter it is masked to force everything outside of -.5 and .5 to zero. If 
    sample_positions is None, then sample_num must be provided. In this case, 
    the filter returned ranges from -.5 to .5 with sample_num of points. 

    The following windows are available: 

    'hann'         von Hann window. Forms a PU. The Hann window has a
                 mainlobe with of 8/N, a PSL of -31.5 dB and decay rate
                 of 18 dB/Octave. 

    'cos'          Cosine window. This is the square root of the Hanning
                 window. The cosine window has a mainlobe width of 6/N,
                 a  PSL of -22.3 dB and decay rate of 12 dB/Octave.
               
    'rec'          Rectangular window. The rectangular window has a
                 mainlobe width of 4/N, a  PSL of -13.3 dB and decay
                 rate of 6 dB/Octave. Forms a PU. Alias: 'square' 

    'tri'          Triangular window.  

    'hamming'      Hamming window. Forms a PU that sums to 1.08 instead
                 of 1.0 as usual. The Hamming window has a
                 mainlobe width of 8/N, a  PSL of -42.7 dB and decay
                 rate of 6 dB/Octave. 

    'blackman'     Blackman window. The Blackman window has a
                 mainlobe width of 12/N, a PSL of -58.1 dB and decay
                 rate of 18 dB/Octave. 

    'blackharr'    Blackman-Harris window. The Blackman-Harris window has 
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

    See also:  gen_filterbank, gen_inv_filterbank
    '''

    # Input argument checking -- either sample_positions or sample_num must be provided
    if sample_positions is not None:
        pass
    elif sample_num != None:
        step = 1/sample_num
        # Two cases if the num of samples is even or odd
        if sample_num % 2 == 0:
            # For even N the sampling interval is [-.5,.5-1/N]
            first_half = np.linspace(0,.5-step,int(sample_num*.5))
            second_half = np.linspace(-.5,-step,int(sample_num*.5))
            sample_positions = np.concatenate((first_half,second_half))
        else: 
            # For odd N the sampling interval is [-.5+1/(2N),.5-1/(2N)]
            first_half = np.linspace(0,.5-.5*step,int(sample_num*.5)+1)
            second_half = np.linspace(-.5+.5*step,-step,int(sample_num*.5))
            sample_positions = np.concatenate((first_half,second_half))
    else:
        print("Error: invalid arguements to window function generator.")
        return None


    # Switch case for possible window names, forcing everything to lowercase
    window_name = window_name.lower()
    if window_name in ['hann','nuttall10']:
        g = .5 + .5*np.cos(2*np.pi*sample_positions)
        
    elif window_name in ['cosine','cos','sqrthann']:
        g = np.cos(np.pi*sample_positions)
        
    elif window_name in ['hamming','nuttall01']:
        g = .54 + .46*np.cos(2*np.pi*sample_positions)
        
    elif window_name in ['square','rec','boxcar']:
        g = np.asarray([int(abs(i) < .5) for i in sample_positions], dtype=np.float64)
        
    elif window_name in ['tri','triangular','bartlett']:
        g = 1-2*abs(sample_positions)
        
    elif window_name in ['blackman']:
        g = .42 + .5*np.cos(2*np.pi*sample_positions) + .08*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['blackharr']:
        g = .35875 + .48829*np.cos(2*np.pi*sample_positions) + .14128*np.cos(4*np.pi*sample_positions) + \
            .01168*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['modblackharr']:
        g = .35872 + .48832*np.cos(2*np.pi*sample_positions) + .14128*np.cos(4*np.pi*sample_positions) + \
            .01168*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall','nuttall12']:
        g = .355768 + .487396*np.cos(2*np.pi*sample_positions) + .144232*np.cos(4*np.pi*sample_positions) + \
            .012604*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall20']:
        g = 3/8 + 4/8*np.cos(2*np.pi*sample_positions) + 1/8*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['nuttall11']:
        g = .40897 + .5*np.cos(2*np.pi*sample_positions) + .09103*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['nuttall02']:
        g = .4243801 + .4973406*np.cos(2*np.pi*sample_positions) + .0782793*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['nuttall30']:
        g = 10/32 + 15/32*np.cos(2*np.pi*sample_positions) + 6/32*np.cos(4*np.pi*sample_positions) + \
            1/32*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall21']:
        g = .338946 + .481973*np.cos(2*np.pi*sample_positions) + .161054*np.cos(4*np.pi*sample_positions) + \
            .018027*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall03']:
        g = .3635819 + .4891775*np.cos(2*np.pi*sample_positions) + .1365995*np.cos(4*np.pi*sample_positions) + \
            .0106411*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['gauss','truncgauss']:
        g = np.exp(-18*sample_positions**2)
        
    elif window_name in ['wp2inp']:
        g = np.exp(np.exp(-2*sample_positions)*25.*(1+2*sample_positions))
        g = g/np.max(g)
    else:
        print("Error: unknown window function name ", window_name, " provided.")
        return None

    # List comprehension makes a mask to force values outside -.5 and .5 to zero
    mask = [int(abs(i) < .5) for i in sample_positions]
    g *= mask
    return g