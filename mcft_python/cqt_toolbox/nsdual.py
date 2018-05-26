from __future__ import print_function,division

import numpy as np


def nsdual(g,shift,M=None):
	'''
	NSDUAL  Canonical dual NSG frame (for painless systems)
	   Usage: gd = nsdual(g,shift,M)
	
	   Input parameters:
	         g         : Cell array of window functions/filters
	         shift     : Vector of time/frequency shifts
	         M         : Number of frequency channels (vector/scalar)
	   Output parameters:
	         gd        : Dual window functions 
	
	   Given a non-stationary Gabor frame specified by the windows g, shift 
	   parameters shift, and channel numbers M, NSDUAL computes the
	   canonical dual frame windows/filters gd by inverting the diagonal of 
	   the frame operator and applying the inverse to g. More explicitly,
	
	      gd{n} = g{n} / ( sum M(l) |g{l}|^2 ), 
	                        l  
	
	   If g, shift, M specify a painless frame, i.e. 
	   SUPP(G{N})  <= M(n)~forall~n and 
	
	      A <= sum ( M(n) |g{n}|^2 ) <= B, for some 0 < A <= B < infty
	            n  
	
	   the computation will result in the canonical dual frame. If  g, 
	   shift, M specify a frame, but the first condition is violated, the 
	   result can be interpreted as a first approximation of the corresponding 
	   canonical dual frame.
	 
	   Note, the time shifts corresponding to the dual window sequence is the
	   same as the original shift sequence and as such already given.
	
	   If g, shift, M is a painless frame, the output can be used for 
	   perfect reconstruction of a signal using the inverse nonstationary 
	   Gabor transform NSIGT.
	 
	   See also:  nsgt, nsigt, nsgt_real, nsigt_real, nsgtf, nsigtf
	'''
	# Argument handling
	if M is None:
		M = np.zeros(len(shift))
		for i in range(len(shift)):
			M[i] = len(g[i])

	if M.size == 1:
		M = M[0]*np.ones(len(shift))

	N = len(shift)
	posit = np.cumsum(shift)
	Ls = posit[N-1]
	posit -= shift[0]
	diag = np.zeros(int(Ls))
	win_range = []

	for i in range(N):
		Lg = len(g[i])
		win_range.append(((posit[i] + np.arange(-1*np.floor(Lg/2),np.ceil(Lg/2))) % Ls).astype(np.int32))
		diag[win_range[i]] += (np.fft.fftshift(g[i])**2)*M[i]

	gd = g
	for i in range(N):
		gd[i] = np.fft.ifftshift(np.fft.fftshift(gd[i])/diag[win_range[i]])
	
	return gd