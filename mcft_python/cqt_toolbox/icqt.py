from __future__ import print_function,division

import numpy as np

from gen_inv_filterbank import gen_inv_filterbank
from apply_inv_filterbank import  apply_inv_filterbank


def icqt(Xcq):
	'''
	ICQT  inverse Constant-Q/Variable-Q transform
	   Usage:  [x gd] = icqt(Xcq, varargin)
	
	   Input parameters:
	         Xcq       : struct obtained by cqt(...)
	         varargin  : Optional input pairs (see table below)
	
	   Output parameters: 
	         x         : reconstructed time domain signal
	         gd        : synthesis filterbank
	
	   See also:  cqt, nsigtf_real, winfuns
	'''
	Xcq['gd'] = gen_inv_filterbank(Xcq['g'],Xcq['shift'],Xcq['M'])

	# We currently assume rasterize is always full
	c = [x for x in Xcq['c']]
	c.insert(0,Xcq['cDC'])
	c.append(Xcq['cNyq'])
	
	x = apply_inv_filterbank(c,Xcq['gd'],Xcq['shift'],Xcq['xlen'],Xcq['phasemode'])
	
	gd = Xcq['gd']
	
	return x, gd