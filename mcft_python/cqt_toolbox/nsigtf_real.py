from __future__ import print_function,division

import numpy as np


def nsigtf_real(c,g,shift,Ls,phasemode):
	'''
	 NSIGTF_REAL  Nonstationary Gabor filterbank synthesis for real signals
	   Usage: fr = nsigtf_real(c,g,shift,M,Ls)
	
	   Input parameters: 
	         c         : Cell array of nonstationary Gabor coefficients
	         g         : Cell array of synthesis filters
	         shift     : Vector of time shifts
	         M         : Number of time channels (vector/scalar)
	         Ls        : Length of the analyzed signal
	         phasemode  :  can be set to (default is 'global')
	                      - 'local':     Zero-centered filtered used
	                      - 'global':    Mapping function used (see reference)
	   Output parameters:
	         fr        : Synthesized real-valued signal (Channels are stored 
	                     in the columns)
	
	   Given the cell array c of nonstationary Gabor filterbank 
	   coefficients, a set of filters g and frequency shifts shift, this 
	   function computes the corresponding nonstationary Gabor filterbank
	   synthesis for real valued signals. 
	
	   Note that, due to the structure of the coefficient array in the real
	   valued setting, all entries g{n} with N > length(c) will be ignored
	   and assumed to be fully supported on the negative frequencies.
	
	   Let P(n)=sum_{l=1}^{n} shift(l), then the synthesis formula reads:
	
	                    N-1 
	       fr_temp(l) = sum sum c{n}(m)g{n}[l-P(n)]*exp(-2*pi*i*(l-P(n))*m/M(n)),
	                    n=0  m
	   
	   for l=0,cdots,Ls-1.  In practice, the synthesis formula is realized 
	   by fft and overlap-add. To synthesize the negative frequencies, 
	   fr_temp is truncated to length floor( Ls/2 )+1. Afterwards 
	   ifftreal implicitly computes the hermite symmetric extension and 
	   computes the inverse Fourier transform, i.e. fr = ifftreal(fr_temp).
	 
	   If a nonstationary Gabor frame was used to produce the coefficients 
	   and g is a corresponding dual frame, this function should perfectly 
	   reconstruct the originally analyzed signal to numerical precision.
	   
	   Multichannel output will save each channel in a column of fr.
	
	   See also:  nsdual, nstight
	'''
	if len(c[0].shape) > 1:
		CH = c[0].shape[1]
	else:
		CH = 1
	N = len(c)	

	posit = np.cumsum(shift)
	NN = posit[-1]
	posit -= shift[0]

	fr = np.zeros((int(NN),int(CH)),dtype=np.complex128)

	for i in range(N):
		Lg = len(g[i])


		win_range = ((posit[i] + np.arange(-1*np.floor(Lg/2),np.ceil(Lg/2))) % NN).astype(np.int32)

		temp = np.fft.fft(c[i],axis=0)*len(c[i])
		
		if phasemode == 'global':
			fsNewBins = len(c[i])
			fkBins = posit[i]
			displace = fkBins - np.floor(fkBins/fsNewBins) * fsNewBins
			temp = np.roll(temp, int(displace))

		first_half = np.arange(len(temp)-np.floor(Lg/2),len(temp))
		second_half = np.arange(np.ceil(Lg/2))
		idx1 = np.concatenate((first_half,second_half))
		temp = temp[np.mod(idx1,len(temp)).astype(np.int32)]
		idx2 = np.concatenate((np.arange(Lg-np.floor(Lg/2),Lg),np.arange(np.ceil(Lg/2)))).astype(np.int32)
		fr[win_range,:] += (temp * g[i][idx2]).reshape(len(temp),1)

	nyqBin = int(np.floor(Ls/2))
	fr[nyqBin+1:] = np.conj(fr[nyqBin - (1- (Ls%2)):0:-1])
	fr = np.real(np.fft.ifft(fr,axis=0))

	return fr