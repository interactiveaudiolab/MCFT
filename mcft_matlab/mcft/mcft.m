function [mcft_out,cqt_params_out,H]=mcft(x,cqt_params_in,filt2d_params)

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
% Sch�rkhuber et al. "A Matlab toolbox for efficient perfect 
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
%        scale_ctrs: vector containing filter centers (along scale axis)
%        rate_ctrs: vector containing filter centers (along rate axis)
%        time_const: time constant of the temporal filter
%
% Outputs:
% mcft_out: 4d matrix containing MCFT coefficients
% cqt_params_out: structure array containing all CQT properties 
%                 (required for reconstruction)
% H: 4d matrix containing the scale-rate filter bank
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Input check

if nargin<3
    set_filt_params = 1;
else
    set_filt_params = 0;
    scale_ctrs = filt2d_params.scale_ctrs;
    rate_ctrs = filt2d_params.rate_ctrs;
    beta = filt2d_params.time_const;
end


%% CQT parameters:

fs = cqt_params_in.fs;   % sample rate of the time signal
fmin = cqt_params_in.fmin; % min frequency of the cqt filter bank 
fmax = cqt_params_in.fmax; % max frequency of the cqt filter bank
fres = cqt_params_in.fres; % number of filters per octave
gamma = cqt_params_in.gamma; % linear-Q factor

%% Time-domain signal to CQT

Xcq = cqt(x,fres, fs, fmin, fmax, 'rasterize', 'full','gamma',gamma);
X = Xcq.c;
[Nf,Nt] = size(X); % number of frequency channels and time frames

% observation: this part is invertible for just one
% signal (without any manipulation) but fails in separation
% reason: phase is not linear over harmonics in cqt
% it's actually periodic (over time) with different period
% for each overtone, when unwrap it over time you can
% clearly see lines with different slopes
% needs more work!
%%%% phase normalization method %%%%
% % unwrap and frequency normalize the phase
% freq_vec=Xcq.fbas;
% fmat=repmat(freq_vec,1,Nt);
% Xph=unwrap(angle(X),[],2)./fmat;
%%% the following line is interesting and worth exploring
% X=abs(X).*exp(1j*abs(X).*angle(X)); %Xph);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cqt_params_out = Xcq;
cqt_params_out = rmfield(cqt_params_out,'c');
cqt_params_out.Nf = Nf;
cqt_params_out.Nt = Nt;

 
%% Parameters of spectro-temporal filters:

x_dur = length(x)/fs; % duration of the signal (in sec)

nfft_s = Nf; % min number of fft points along the scale axis 
nfft_r = Nt; % min number of fft points along the rate axis 

samprate_spec=fres; % sampling rate of the spectral filter (in samples per octave)
samprate_temp=floor(Nt/x_dur); % sampling rate of the temporal filter (in frames per sec)

% set filter scales and rates if not provided
if set_filt_params   
   [scale_ctrs,rate_ctrs]=filt_default_centers(nfft_s,nfft_r,samprate_spec,samprate_temp);
    beta=1; % time constant of the temporal filter
end
    
% concatenate filter parameters into one structure array
H_params=struct('samprate_spec',samprate_spec,'samprate_temp',samprate_temp,'time_const',beta);
    

%% Spectro-temporal filter bank

disp('Computing the filterbank...');
%%% the following line is interesting and worth exploring
% aa = abs(X).*exp(1j*abs(X).*angle(X));
[~,H]=gen_fbank_hsr(scale_ctrs,rate_ctrs,nfft_s,nfft_r,H_params,X); 

%% CQT to MCFT

disp('Computing the transform...');
mcft_out=cqt_to_mcft(X,H);

end



