from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt

def gen_filter(window_name, sample_positions=None, window_len=None):
    '''
    WINFUNS  Window function generator  
       Usage:  g = winfuns(name,x)
               g = winfuns(name,N,L)
               g = winfuns(name,N)
       
       Input parameters: 
             name      : String containing the window name
             x         : Vector of sampling positions
             N         : Window support (in samples)
             L         : Output length (in samples)
       Output parameters:
             g         : Output window
       
       This function serves to compute a variety of standard and some more 
       exotic window functions. Most of the functions used are detailed and 
       discussed in classical papers (see references below), but several are
       included for special purposes in the toolbox only.
       
       Given a character string name containing the name of the desired
       window function, the function offers 2 modes of operation. If the 
       second input parameter is a vector x of sampling values, then the
       specified function is evaluated at the given points. If a window length
       N and optionally a signal length L are supplied, a symmetric, 
       whole-point centered window with a support of N samples is produced 
       and, given L, zero-extended to length L.

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

       See also:  nsgcqwin, nsgwvltwin, nsgerbwin
    '''
    if sample_positions is not None:
        pass  # TODO: should probably make a check that input is a numpy array
    elif window_len != None:
        step = 1/window_len
        if window_len % 2 == 0:
            first_half = np.linspace(0,.5-step,int(window_len*.5))
            second_half = np.linspace(-.5,-step,int(window_len*.5))
            sample_positions = np.concatenate((first_half,second_half))
        else: 
            first_half = np.linspace(0,.5-.5*step,int(window_len*.5)+1)
            second_half = np.linspace(-.5+.5*step,-step,int(window_len*.5))
            sample_positions = np.concatenate((first_half,second_half))
    else:
        print("Error: invalid arguements to window function generator.")
        return None

    # TODO: Missing check to make sure x is a column vector:
    # this may or may not be necessary depending on how exactly
    # numpy handles vectors as opposed to MATLAB

    if window_name in ['Hann','hann','nuttall10','Nuttall10']:
        g = .5 + .5*np.cos(2*np.pi*sample_positions)
        
    elif window_name in ['Cosine','cosine','cos','Cos','sqrthann','Sqrthann']:
        g = np.cos(np.pi*sample_positions)
        
    elif window_name in ['hamming','nuttall01','Hamming','Nuttall01']:
        g = .54 + .46*np.cos(2*np.pi*sample_positions)
        
    elif window_name in ['square','rec','Square','Rec','boxcar','Boxcar']:
        g = np.asarray([int(abs(i) < .5) for i in sample_positions], dtype=np.float64)
        
    elif window_name in ['tri','triangular','bartlett','Tri','Triangular','Bartlett']:
        g = 1-2*abs(sample_positions)
        
    elif window_name in ['blackman','Blackman']:
        g = .42 + .5*np.cos(2*np.pi*sample_positions) + .08*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['blackharr','Blackharr']:
        g = .35875 + .48829*np.cos(2*np.pi*sample_positions) + .14128*np.cos(4*np.pi*sample_positions) + \
            .01168*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['modblackharr','Modblackharr']:
        g = .35872 + .48832*np.cos(2*np.pi*sample_positions) + .14128*np.cos(4*np.pi*sample_positions) + \
            .01168*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall','nuttall12','Nuttall','Nuttall12']:
        g = .355768 + .487396*np.cos(2*np.pi*sample_positions) + .144232*np.cos(4*np.pi*sample_positions) + \
            .012604*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall20','Nuttall20']:
        g = 3/8 + 4/8*np.cos(2*np.pi*sample_positions) + 1/8*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['nuttall11','Nuttall11']:
        g = .40897 + .5*np.cos(2*np.pi*sample_positions) + .09103*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['nuttall02','Nuttall02']:
        g = .4243801 + .4973406*np.cos(2*np.pi*sample_positions) + .0782793*np.cos(4*np.pi*sample_positions)
        
    elif window_name in ['nuttall30','Nuttall30']:
        g = 10/32 + 15/32*np.cos(2*np.pi*sample_positions) + 6/32*np.cos(4*np.pi*sample_positions) + \
            1/32*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall21','Nuttall21']:
        g = .338946 + .481973*np.cos(2*np.pi*sample_positions) + .161054*np.cos(4*np.pi*sample_positions) + \
            .018027*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['nuttall03','Nuttall03']:
        g = .3635819 + .4891775*np.cos(2*np.pi*sample_positions) + .1365995*np.cos(4*np.pi*sample_positions) + \
            .0106411*np.cos(6*np.pi*sample_positions)
        
    elif window_name in ['gauss','truncgauss','Gauss','Truncgauss']:
        g = np.exp(-18*sample_positions**2)
        
    elif window_name in ['wp2inp','Wp2inp']:
        g = np.exp(np.exp(-2*sample_positions)*25.*(1+2*sample_positions))
        g = g/np.max(g)
    else:
        print("Error: unknown window function name ", window_name, " provided.")
        return None

    mask = [int(abs(i) < .5) for i in sample_positions]
    g *= mask
    return g