import numpy as np 

def nsgcqwin(fmin,fmax,bins,sr,Ls,min_win=4,bwfac=1,fractional=0,winfun='hann',gamma=0):
	'''
	NSGCQWIN  Constant-Q/Variable-Q dictionary generator
	   Usage:  [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,Ls,varargin)
	           [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,Ls)
	
	   Input parameters:
	         fmin      : Minimum frequency (in Hz)
	         fmax      : Maximum frequency (in Hz)
	         bins      : number of bins per octave
	         sr        : Sampling rate (in Hz)
	         Ls        : Length of signal (in samples)
	         varargin  : Optional input pairs (see table below)
	   Output parameters: 
	         g         : Cell array of constant-Q/variable-Q filters
	         shift     : Vector of shifts between the center frequencies
	         M         : Vector of lengths of the window functions
	
	   Create a nonstationary Gabor filterbank with constant or varying 
	   Q-factor and relevant frequency range from fmin to fmax. To allow
	   for perfect reconstruction, the frequencies outside that range will be
	   captured by 2 additional filters placed on the zero and Nyquist
	   frequencies, respectively.
	
	   The Q-factor (quality factor) is the ratio of center frequency to
	   bandwidth cent_freq/bandwidth.
	
	
	   For more details on the construction of the constant-Q nonstationary 
	   Gabor filterbank, please check the reference.
	   
	   Optional input arguments arguments can be supplied like this:
	
	       nsgcqwin(fmin,fmax,bins,sr,Ls,'min_win',min_win)
	
	   The arguments must be character strings followed by an
	   argument:
	
	     'min_win',min_win  Minimum admissible window length (in samples) 
	
	     'bwfac',bwfac            Channel numbers M are rounded to multiples 
	                              of this
	
	     'fractional',fractional  Allow fractional shifts and bandwidths
	
	     'winfun',winfun          String containing the desired window 
	                              function name
	
	     'gamma':      the bandwidth of each filter is given by
	                            Bk = 1/Q * fk + gamma,
	                   where fk is the filters center frequency, Q is fully
	                   determined by the number of bins per octave and gamma
	                   is a bandwidth offset. If gamma = 0 the obtained
	                   filterbank is constant-Q. Setting gamma > 0 time
	                   resolution towards lower frequencies can be improved
	                   compared to the constant-Q case (e.g. ERB proportional
	                   bandwidths). See reference for more information.
	
	   See also: nsgtf_real, winfuns
	
	   References:
	     C. Schörkhuber, A. Klapuri, N. Holighaus, and M. Dörfler. A Matlab 
	     Toolbox for Efficient Perfect Reconstruction log-f Time-Frequecy 
	     Transforms.
	
	     G. A. Velasco, N. Holighaus, M. DÃ¶rfler, and T. Grill. Constructing an
	     invertible constant-Q transform with non-stationary Gabor frames.
	     Proceedings of DAFX11, Paris, 2011.
	     
	     N. Holighaus, M. DÃ¶rfler, G. Velasco, and T. Grill. A framework for
	     invertible, real-time constant-q transforms. Audio, Speech, and
	     Language Processing, IEEE Transactions on, 21(4):775-785, April 2013.
	     
	
	   Url: http://nsg.sourceforge.net/doc/generators/nsgcqwin.php
	'''
	nyquist = sr/2.
	if fmax > nyquist:
		fmax = nyquist

	fftres = sr/Ls
	b = np.floor(bins * np.log2(fmax/fmin))
	fbas = fmin * 2**(np.asarray(range(0,b+1))/bins) 

	Q = 2**(1/bins) - 2**(-1/bins)
	cqtbw = Q*fbas + gamma

	nonzeroIndices = np.nonzero([int(fbas[i]+cqtbw[i]/2>nyquist) for i in range(fbas.shape[0])])
	if nonzeroIndices[0].size > 0:
		fbas = fbas[:nonzeroIndices[0][0]]
		cqtbw = cqtbw[:nonzeroIndices[0]]

	nonzeroIndices = np.nonzero([int(fbas[i]-cqtbw[i]/2<0) for i in range(fbas.shape[0])])
	if nonzeroIndices[0].size > 0:
		fbas = fbas[nonzeroIndices[0][-1]:]
		cqtbw = cqtbw[nonzeroIndices[0][-1]:]
		print "fmin set to ", None, " Hz!" # Computation 

	Lfbas = fbas.shape[0]
	fbas = np.concatenate(([0],fbas,[nyquist],sr-np.flip(fbas,0)))
	bw = np.concatenate((2*fmin, cqtbw, fbas[Lfbas+2]-fbas[Lfbas], np.flip(cqtbw,0)))
	bw /= fftres
	fbas /= fftres

	posit = np.zeros(fbas.size)
	posit[:Lfbas+2] = np.floor(fbas[:Lfbas+2])
	posit[Lfbas+2:] = np.ceil(fbas[Lfbas+2:])

	shift = np.concatenate((-1*posit[-1] % Lfbas, np.diff(posit)))

	if fractional:
		corr_shift = fbas-posit
		M = np.ceil(bw+1)
	else:
		bw = np.round(bw)
		M = bw

	for i in range(2*(Lfbas+1)):
		if bw[i] < min_win:
			bw[i] = min_win
			M[i] = bw[i]

	if fractional:
		