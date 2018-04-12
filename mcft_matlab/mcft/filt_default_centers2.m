function [S_vec,R_vec]=filt_default_centers2(SRF,FPS)

% This function computes the default set of filter centers
%
% Inputs:
% SRF: sampling rate of the spectral filter (in samples per octave)
% FPS: sampling rate of the temporal filter (in frames per sec)
%
% Outputs:
% SV: vector containing filter centers (along scale axis)
% RV: vector containing filter centers (along rate axis)
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% set filter scales
    log2_slow=-3; % center of the low-pass filter
    log2_sbandmax=nextpow2(SRF/2)-1; % center of the highest band-pass filter
    log2_shigh=log2((2^log2_sbandmax + SRF/2)/2); % center of the high-pass filter
    S_vec=2.^[(log2_slow:1:log2_sbandmax),log2_shigh]; 
    
%% set filter rates
    log2_rlow=-3; % center of the low-pass filter
    log2_rbandmax=nextpow2(FPS/2)-1; % center of the highest band-pass filter
    log2_rhigh=log2((2^log2_rbandmax + FPS/2)/2); % center of the high-pass filter
    R_vec=2.^[(log2_rlow:0.5:log2_rbandmax),log2_rhigh];
end