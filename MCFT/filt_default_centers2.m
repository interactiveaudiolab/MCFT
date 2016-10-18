function [S_vec,R_vec]=filt_default_centers2(nfft_s,nfft_r,SRF,FPS)

% This function computes the default set of filter centers
%
% Inputs:
% nfft_s: number of fft points on the scale axis
% nfft_r: number of fft points on the rate axis
% SRF: sampling rate of the spectral filter (in samples per octave)
% FPS: sampling rate of the temporal filter (in frames per sec)
%
% Outputs:
% SV: vector containing filter centers (along scale axis)
% RV: vector containing filter centers (along rate axis)
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% set filter scales
    smin=SRF/nfft_s; % center of the low-pass filter
    log2_sbandmax=ceil(log2(nfft_s/2)-1); % center of the highest band-pass filter
    shigh=(smin*2^log2_sbandmax + SRF/2)/2; 
    shigh=smin*floor(shigh/smin);   % center of the high-pass filter
    S_vec=[smin*2.^(-1:log2_sbandmax),shigh];
    
%% set filter rates
    rmin=FPS/nfft_r; % center of the low-pass filter
    log2_rbandmax=ceil(log2(nfft_r/2)-1); % center of the highest band-pass filter
    rhigh=(rmin*2^log2_rbandmax + FPS/2)/2; 
    rhigh=rmin*floor(rhigh/rmin);   % center of the high-pass filter
    R_vec=[rmin*2.^(-1:0.5:log2_rbandmax),rhigh];
 
end