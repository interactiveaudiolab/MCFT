from __future__ import print_function, division

import numpy as np

from nsgcqwin import nsgcqwin
from nsgtf_real import  nsgtf_real


def cqt(x, B, fs, fmin, fmax,
			rasterize='full', phasemode='global', format='sparse',
			gamma=0, normalize='sine', win_fun='hann'):
	'''
	%CQT  Constant-Q/Variable-Q transform
	%   Usage:  Xcq = cqt(x, B, fs, fmin, fmax, varargin)
	%
	%   Input parameters:
	%         x         : input signal
	%         B         : number of bins per octave
	%         fs        : sampling frequency
	%         fmin      : lowest frequency to be analyzed
	%         fmax      : highest frequency to be analyzed
	%         varargin  : Optional input pairs (see table below)
	%
	%   Output parameters: 
	%         Xcq       : Struct consisting of 
	%           .c           : CQT coefficients
	%           .cDC         : transform coefficients for f = 0
	%           .cNyq        : transform coefficients for fs/2
	%           .g           : cell array of analysis filters
	%           .shift       : center frequencies of analysis filters
	%           .M           : bandwidth of analysis filters
	%           .xlen        : length of input signal
	%           .phasemode   : 'local'  -> zero-centered filtered used
	%                        : 'global' -> mapping function used
	%           .rast        : time-frequency plane sampling scheme (full,
	%                          piecewise, none)
	%           .fmin
	%           .fmax
	%           .B       
	%           .format      : eighter 'cell' or 'matrix' (only applies for
	%                          piecewise rasterization)
	%   
	%   Optional input arguments arguments can be supplied like this:
	%
	%       Xcq = cqt(x, B, fs, fmin, fmax, 'rasterize', 'piecewise')
	%
	%   The arguments must be character strings followed by an
	%   argument:
	%
	%     'rasterize':  can be set to (default is 'full');
	%           - 'none':      Hop sizes are distinct for each frequency
	%                          channel. Transform coefficients will be
	%                          presented in a cell array.
	%           - 'full':      The hop sizes for all freqency channels are 
	%                          set to the smallest hop size in the representa-
	%                          tion. Transform coefficients will be presented 
	%                          in matrix format.
	%           - 'piecewise': Hop sizes will be rounded down to be a power-of-
	%                          two integer multiple of the smallest hop size in
	%                          the representation. Coefficients will be 
	%                          presented either in a sparse matrix or as cell 
	%                          arrays (see 'format' option)
	%
	%     'phasemode':  can be set to (default is 'global')
	%           - 'local':     Zero-centered filtered used
	%           - 'global':    Mapping function used (see reference)
	%
	%     'format':     applies only for piecewise rasterization               
	%           - 'sparse':   Coefficients will be presented in a sparse matrix 
	%           - 'cell':     Coefficients will be presented in a cell array
	%
	%     'gamma':      the bandwidth of each filter is given by
	%                            Bk = 1/Q * fk + gamma,
	%                   where fk is the filters center frequency, Q is fully
	%                   determined by the number of bins per octave and gamma
	%                   is a bandwidth offset. If gamma = 0 the obtained
	%                   filterbank is constant-Q. Setting gamma > 0 time
	%                   resolution towards lower frequencies can be improved
	%                   compared to the constant-Q case (e.g. ERB proportional
	%                   bandwidths). See reference for more information.
	%     'normalize':  coefficient normalization
	%          - 'sine':    Filters are scaled such that a sinusoid with
	%                       amplitude A in time domain will exhibit the same
	%                       amplitude in the time-frequency representation.
	%          - 'impulse': Filters are scaled such that an impulse in time
	%                       domain will exhibit a flat response in the
	%                       time-frequency representation (in the frame that 
	%                       centers the impulse)
	%          - 'none':      ...
	%     'winfun':        defines the window function that is used for filter
	%                   design. See winfuns for more information.
	%
	%   See also:  nsgtf_real, winfuns
	'''

	g,shift,M = nsgcqwin(fmin,fmax,B,fs,len(x), winfun=win_fun, gamma=gamma)

	total_bins = int(len(M)/2 -1)
	fbas = fs * np.cumsum(shift[1:]) / len(x)
	fbas = fbas[:total_bins+1]

	# Assumes rasterize is full always
	M[1:total_bins+2] = M[total_bins+1]
	M[total_bins+2:] = M[total_bins:0:-1]

	if normalize in ['sine','Sine','SINE','sin']:
		normFacVec = 2 * M[:total_bins+3]/len(x)
	elif normalize in ['impulse','Impulse','IMPULSE','imp']:
		win_lengths = np.zeros(len(g))
		for i in range(len(g)):
			win_lengths[i] = len(g[i])
		normFacVec = 2 * M[:total_bins+3]/win_lengths
	elif normalize in ['none','None','NONE','no']:
		normFacVec = np.ones(total_bins+2)
	else:
		print("Unknown normalization method")

	normFacVec = np.concatenate((normFacVec,normFacVec[len(normFacVec)-2:0:-1]))
	import pdb; pdb.set_trace()

	for i in range(len(normFacVec)):
		g[i] *= normFacVec[i]

	c = nsgtf_real(x,g,shift,phasemode,M)

	# Assume rasterize is full always
	cDC = c[0]
	cNyq = c[total_bins+1]
	c = c[1:total_bins+1]

	results = {'c':c,'g':g,'shift':shift,'M':M,'xlen':len(x),
		'phasemode':phasemode,'rast':rasterize,'fmin':fmin,'fmax':fmax,
		'B':B,'cDC':cDC,'cNyq':cNyq,'format':outputFormat,'fbas':fbas}

	return results