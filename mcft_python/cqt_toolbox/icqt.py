from __future__ import print_function,division

import numpy as np

from nsdual import nsdual
from nsigtf_real import  nsigtf_real

from scipy.io import loadmat
import matplotlib.pyplot as plt

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
	Xcq['gd'] = nsdual(Xcq['g'],Xcq['shift'],Xcq['M'])

	# We currently assume rasterize is always full
	c = [x for x in Xcq['c']]
	c.insert(0,Xcq['cDC'])
	c.append(Xcq['cNyq'])

	matlab = loadmat('inputs.mat')
	import pdb; pdb.set_trace()
	
	x = nsigtf_real(c,Xcq['gd'],Xcq['shift'],Xcq['xlen'],Xcq['phasemode'])
	
	gd = Xcq['gd'][0]

	return x, gd