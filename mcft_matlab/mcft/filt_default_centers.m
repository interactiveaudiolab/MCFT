function [scale_ctrs,rate_ctrs] = filt_default_centers(nfft_s,nfft_r,samprate_spec,samprate_temp)

% This function computes the default set of filter centers
%
% Inputs:
% nfft_s: number of fft points on the scale axis
% nfft_r: number of fft points on the rate axis
% samprate_spec: sampling rate of the spectral filter (in cycles per octave)
% samprate_temp: sampling rate of the temporal filter (in cycles per second)
%
% Outputs:
% scale_ctrs: vector containing filter centers (along scale axis)
% rate_ctrs: vector containing filter centers (along rate axis)
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% set filter scales
    s_res = samprate_spec / nfft_s; % resolution of the s grid
    cs_low = s_res / 2; % center of the low-pass filter
    log2_cs_band_max = nextpow2(samprate_spec / 2) - 1; % center of the highest band-pass filter
    cs_band = 2.^(0:1:log2_cs_band_max); % centers of band-pass filters
    cs_high = (cs_band(end) + 3*samprate_spec / 2) / 4; % center of the high-pass filter
    
    % move the centers to the nearest frequencies of analysis
    cs_band = round(cs_band/s_res)*s_res;
    cs_high = round(cs_high/s_res)*s_res;
    
    % concatenate all center into one vector
    scale_ctrs = [cs_low,cs_band,cs_high];
    
    % remove repeated values (sometimes happens due to rounding)
    scale_ctrs=unique(scale_ctrs);
   
%% set filter rates
    r_res = samprate_temp / nfft_r; % resolution of the r grid
    cr_low = r_res / 2; % center of the low-pass filter
    log2_cr_band_max = nextpow2(samprate_temp / 2) - 1; % center of the highest band-pass filter
    cr_band = 2.^(0:1:log2_cr_band_max); % centers of band-pass filters
    cr_high = (cr_band(end) + 1*samprate_temp / 2) / 2; % center of the high-pass filter
    
    % move the centers to the nearest frequencies of analysis
    cr_band = round(cr_band/r_res)*r_res;
    cr_high = round(cr_high/r_res)*r_res;
    
    % concatenate all center into one vector
    rate_ctrs = [cr_low,cr_band,cr_high];
    
    % remove repeated values (sometimes happens due to rounding)
    rate_ctrs=unique(rate_ctrs);
 
end