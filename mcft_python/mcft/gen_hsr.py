import numpy as np
from collections import namedtuple

def gen_hsr(scale_ctr, rate_ctr, scale_params, rate_params, filt_dir):
    """
    This function generates a 2D-impulse response in the time-frequency
    domain with dilation factors S and R: h(omega, tau; S, R)

    Inputs:
    scale_ctr: filter center along the scale axis
    rate_ctr: filter center along the rate axis  
    scale_params: named tuple containing the parameters of the spectral filter
             hslen: length of the spectral filter impulse response
             samprate_spec: sample rate along the freq. axis (cyc/oct)
             type: string argument indicating the filter type
                  ('bandpass','lowpass','highpass')

    Example: namedtuple(scale_params,['length','sample_rate','type'])

    rate_params: named tuple containing the parameters of the temporal filter
             time_const: exponent coefficient of the exponential term
             hrlen: length of the temporal filter impulse response
             samprate_temp: sample rate along the time axis
             type: string argument indicating the type of the filter
                   ('bandpass','lowpass','highpass')

    Example: namedtuple(rate_params,['length','sample_rate','type','time_const'])

    filt_type: string type, determines the moving direction of the filter 
               'none' (full s-r domain)
               'up' (upward analytic, nonzero over 2nd and 3rd quadrants)
               'down' (downward analytic, nonzero over 1st and 4th quadrants)
    
    Outputs
    hsr: impulse response of the 2D filter (in time-frequency domain)
    Hsr: 2D-Fourier transform of the filter impulse response 
    
    Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)
            Trent Cwiok (cwiok@u.northwestern.edu)
    """
    ### Frequency and time vectors
    # Zero-pad filters to next even number
    hs_len = scale_params.length + (scale_params.length % 2)
    hr_len = rate_params.length + (rate_params.length % 2)

    # Frequency and time vectors
    ##### Ask Fatemeh -- is the -1 still correct even with Matlab 1-indexing
    freq_vec = np.arange(scale_params.length)/scale_params.sample_rate
    time_vec = np.arange(scale_params.length)/scale_params.sample_rate

    ### Impulse response of the original scale filter: Gaussian
    hs = scale_ctr*(1-2*np.sqaure(scale_ctr*np.pi*freq_vec))
            * np.exp(-(np.sqaure(scale_ctr*freq_vec*np.pi)))
    hs = np.append(hs[:(hs_len/2)+1],hs[hs_len/2:0:-1])

    ### Impulse response of the original rate filter
    hr = rate_ctr*np.sqaure(rate_ctr*time_vec)
            * np.exp(-1*time_vec*rate_params.time_const*rate_ctr)
            * np.sin(2*np.pi*rate_ctr*time_vec)

    ### Scale Response (Fourier transform of the scale impulse response)
    # band-pass scale filter
    Hs = np.absolute(np.fft.fft(hs,hs_len))
    # band-pass rate filter
    Hr = np.fft.fft(hr, hr_len)

    ### Low/high pass filtering
    if scale_params.type != 'bandpass':
    	Hs1 = Hs[:hs_len/2+1]
    	Hs1 /= np.max(Hs1)
    	Hs2 = Hs[hs_len/2+1:]
    	Hs2 /= np.max(Hs2)
    	max_idx1 = np.argmax(Hs1), max_idx2 = np.argmax(Hs2)

    	if scale_params.type == 'lowpass':
    		Hs1[:max_idx1] = 1
    		Hs2[max_idx2+1:] = 1
    	elif scale_params.type == 'highpass':
    		Hs1[max_idx1+1:] = 1
    		Hs2[:max_idx2] = 1

    	Hs = np.append(Hs1, Hs2)

    if rate_params.type != 'bandpass':
    	Hr_ph = np.angle(Hr)
    	Hr_mag = np.absolute(Hr)

    	Hr_mag1 = Hr_mag[:hr_len/2+1]
    	Hr_mag1 /= np.max(Hr_mag1)
    	Hr_mag2 = Hr_mag[hr_len/2+1:]
    	Hr_mag2 /= np.max(Hr_mag2)
    	max_idx1 = np.argmax(Hr_mag1), max_idx2 = np.argmax(Hr_mag2)

    	if rate_params.type == 'lowpass':
    		Hr_mag1[:max_idx1] = 1
    		Hr_mag2[max_idx2+1:] = 1
    	elif rate_params.type == 'highpass':
    		Hr_mag1[max_idx1+1:] = 1
    		Hr_mag2[:max_idx2] = 1

    	Hr_mag = np.append(Hr_mag1, Hr_mag2)
    	Hr = Hr_mag * np.exp(1j*Hr_ph)

    Hsr_full = np.multiply(Hs, Hr.T)  ## Not sure if this is actually correct
    Hsr_full_mag = abs(Hsr_full)
    Hsr_full_mag /= np.max(Hsr_full_mag)
    Hsr_full = Hsr_full_mag*np.exp(1j*np.angle(Hsr_full))

    if filt_dir == 'up':
    	Hsr_up = Hasr_full
    	Hsr_up[:hs_len/2+1,:hr_len/2+1] = 0
    	Hsr_up[hs_len/2+1:,hr_len/2+1:] = 0
    	Hsr = Hsr_up

    elif filt_dir == 'down':
    	Hsr_down = Hsr_full
    	Hsr_down[:hs_len/2+1,hr_len/2+1:] = 0
    	Hsr_down[hs_len/2+1:,:hr_len/2+1] = 0
    	Hsr = Hsr_down

    else:
    	Hsr = Hsr_full

    hsr = np.fft.ifft2(Hsr)
    if max(imag(hsr))<1e-8:
    	hsr=real(hsr)

    return hsr, Hsr