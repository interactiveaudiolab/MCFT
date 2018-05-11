import numpy as np 


def nsgtf_real(f,g,shift,M=None,phasemode):
	'''
	%NSGTF_REAL  Nonstationary Gabor filterbank for real signals
	%   Usage: [c,Ls] = nsgtf_real(f,g,shift,M, phasemode)
	%          c = nsgtf_real(...)
	%
	%   Input parameters: 
	%         f         : A real-valued signal to be analyzed (For multichannel
	%                     signals, input should be a matrix which each
	%                     column storing a channel of the signal).
	%         g         : Cell array of analysis filters
	%         shift     : Vector of frequency shifts
	%         M         : Number of time channels (optional).
	%                     If M is constant, the output is converted to a
	%                     matrix   ##TC: MUST BE A NUMPY NDARRAY
	%         phasemode : 'local': zero-centered filtered used
	%                     'global': mapping function used (see cqt)
	%   Output parameters:
	%         c         : Transform coefficients (matrix or cell array)
	%         Ls        : Original signal length (in samples)
	%
	%   Given the cell array g of windows, the time shift vector shift, and
	%   channel numbers M, NSGTF_REAL computes the corresponding 
	%   nonstationary Gabor filterbank of f, using only the filters with at 
	%   least partially supported on the positive frequencies. Let 
	%   P(n)=sum_{l=1}^{n} shift(l), then the output 
	%   c = NSGTF_REAL(f,g,shift,M) is a cell array with 
	%
	%              Ls-1                                      
	%      c{n}(m)= sum fft(f)(l)*conj(g{n}(l-P(n)))*exp(2*pi*i*(l-P(n))*m/M(n))
	%               l=0                                      
	%
	%   where m runs from 0 to M(n)-1 and n from 1 to N, where
	%   g{N} is the final filter at least partially supported on the
	%   positive frequencies. All filters in g, shift that are completely
	%   supported on the negative frequencies are ignored.
	%
	%   For more details, see NSGTF.
	%
	%   See also:  nsigtf_real, nsdual, nstight
	'''
	N = len(shift)
	if M is None:
		M = np.zeros(N)
		for i in range(N):
			M[i] = len(g[i])

	if M.size = 1:
		M = M[0]*np.ones(N)

	f = np.fft.fft(f)

	posit = np.cumsum(shift)-shift[0]

	fill = sum(shift)-Ls
	padding = np.zeros((fill,CH))
	f = np.vstack((f,padding))

	Lg = np.zeros(len(g))
	for i in range(len(g)):
		Lg[i] = len(g[i])

	N = [int((posit-np.floor(Lg/2)) <= (Ls+fill)/2)]
	nonzeroIndex = np.nonzero(N)[0][-1]
	c = []

	for i in range(nonzeroIndex):
		idx = np.concatenate((range(int(np.ceil(Lg[i]/2)),int(Lg[i])),range(int(np.ceil(Lg[i]/2)))))
		