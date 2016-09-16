function [Z,cqt_params_out,H]=mcft(x,cqt_params_in,filt2d_params)

% This function receives a time domian audio signal and returns its 
% Multi-resolution Common Fate Transform (MCFT). 
% The intermediary time-frequency domain representation is the 
% Constant-Q Transform (CQT) of the audio signal, which is computed 
% using the invertible and optimized CQT implementation proposed by 
% Schorkhuber et al.:
%
% Toolbox webpage: 
% http://www.cs.tut.fi/sgn/arg/CQT/
% Reference: 
% Schörkhuber et al. "A Matlab toolbox for efficient perfect 
% reconstruction time-frequency transforms with log-frequency resolution."
% Audio Engineering Society Conference: 53rd International Conference: 
% Semantic Audio. Audio Engineering Society, 2014.
%
% Inputs:
% x: vector containing samples of the time-domain signal
% cqt_params_in: structure array containing CQT parameters:
%        fs: sampling rate of the audio signal
%        fmin: minimum frequency of analysis
%        fmax: maximum frequency of analysis
%        fres: frequency resolution (# of bins per octave)
% filt2d_params (optional): structure array containing parameters of 
%        2d (spectro-temporal) filters: 
%        S_vec: vector containing filter centers (along scale axis)
%        R_vec: vector containing filter centers (along rate axis)
%        time_const: time constant of the temporal filter
%
% Outputs:
% Z: 4d matrix containing MCFT coefficients
% cqt_params_out: structure array containing all CQT properties 
%                 (required for reconstruction)
% H: 4d matrix containing the scale-rate filter bank
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Input check

if nargin<3
    set_filt_params=1;
else
    set_filt_params=0;
    SV=filt2d_params.S_vec;
    RV=filt2d_params.R_vec;
    beta=filt2d_params.time_const;
end

%% CQT parameters:

fs=cqt_params_in.fs;
fmin=cqt_params_in.fmin;
fmax=cqt_params_in.fmax;
fres=cqt_params_in.fres;

%% Time-domain signal to CQT

Xcq = cqt(x,fres, fs, fmin, fmax, 'rasterize', 'full','gamma',0);
X=Xcq.c;

[Nf,Nt]=size(X); % number of frequency channels and time frames
cqt_params_out=Xcq;
cqt_params_out=rmfield(cqt_params_out,'c');
cqt_params_out.Nf=Nf;
cqt_params_out.Nt=Nt;

%% Parameters of spectro-temporal filters:

x_dur=length(x)/fs; % duration of the signal (in sec)

nfft_s=Nf; % number of fft points along the scale axis
nfft_r=Nt; % number of fft points along the rate axis 

SRF=fres; % sampling rate of the spectral filter (in samples per octave)
FPS=floor(Nt/x_dur); % sampling rate of the temporal filter (in frames per sec)

% set filter scales and rates if not provided
if set_filt_params   
    % set filter scales
    log2_slow=-3; % center of the low-pass filter
    log2_sbandmax=nextpow2(SRF/2)-1; % center of the highest band-pass filter
    log2_shigh=log2(2^log2_sbandmax + (SRF/2 - 2^log2_sbandmax)/2); % center of the high-pass filter
    SV=2.^[(log2_slow:1:log2_sbandmax),log2_shigh]; 
    
    % set filter rates
    log2_rlow=-3; % center of the low-pass filter
    log2_rbandmax=nextpow2(FPS/2)-1; % center of the highest band-pass filter
    log2_rhigh=log2(2^log2_rbandmax + (FPS/2 - 2^log2_rbandmax)/2); % center of the high-pass filter
    RV=2.^[(log2_rlow:0.5:log2_rbandmax),log2_rhigh];
    
    beta=1; % time constant of the temporal filter
end
    
% concatenate filter parameters into one structure array
H_params=struct('ripple_freq',SRF,'frame_per_sec',FPS,'time_const',beta);
    
%% Spectro-temporal filter bank

disp('Computing the filterbank...');
H=gen_fbank_hsr(SV,RV,nfft_s,nfft_r,H_params,X); 

%% CQT to MCFT

disp('Computing the transform...');
Z=cqt_to_mcft(X,H);

end



