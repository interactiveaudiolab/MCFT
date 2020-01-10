function [mcft_out,cqt_params_out,fbank_sr_domain] = mcft_refactored(signal,cqt_params_in,varargin)

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
% signal: vector containing samples of the time-domain signal
% cqt_params_in: structure array containing CQT parameters:
%        fs: sampling rate of the audio signal
%        fmin: minimum frequency of analysis
%        fmax: maximum frequency of analysis
%        fres: frequency resolution (# of bins per octave)
%        gamma: linear-Q factor
%
% Optional input arguments can be supplied like this:
%      mcft(signal,cqt_params_in,'del_cqt_phase',del_cqt_phase)
% 
% The arguments must be character strings followed by the argument value:
% 'del_cqt_phase',del_cqt_phase: boolean, indicating whether the cqt phase 
%                 is to be deleted or included in the 2d filtering process. 
%                 If the phase is included, the filterbank will be modulated
%                 with the CQT phase.
%                 If the phase is deleted the origianl filterbank will be 
%                 applied to the magnitude CQT.
%                 Note: the output representation is invertible only if the
%                       phase is included (default = 0)
% 
% 'predef_scale_ctrs',predef_scale_ctrs: vector containing filter centers 
%                                       (along scale axis)
% 'predef_rate_ctrs',predef_rate_ctrs: vector containing filter centers 
%                                      (along rate axis)
% 'time_const',time_const: time constant of the temporal filter
%                          default = 1
% 'scale_res',scale_res: number of bins per octvave on the scale axis
%                        default = 1
% 'scale_max',scale_max: the center of the highest bandpass filter
%                        default: 2^(nextpow2(samprate_spec/2)-1) 
%                        (last power2 number before the Nyquist rate)
%                        e.g., scale_max = 8 if samp rate = 24 
% 'scale_min',scale_min: the center of the lowest bandpass filter
%                        default: 2^0                         
% 'rate_res',rate_res: number of bins per octave on the rate axis
%                      default = 1  
% 'rate_max',rate_max: the center of the highest pandpass filter
%                        default: 2^(nextpow2(samprate_temp/2)-1) 
%                        (last power2 number before the Nyquist rate)
%                        e.g., rate_max = 64 if samprate = 200 
% 'rate_min',rate_min: the center of the lowest filter
%                        default: 2^0
%
% Outputs:
% mcft_out: 4d matrix containing MCFT coefficients
% cqt_params_out: structure array containing all CQT properties 
%                 (required for reconstruction)
% fbank_sr_domain: 4d matrix containing the scale-rate filter bank
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Input check

% extract optional inputs
func_inputs = inputParser;
addParameter(func_inputs,'del_cqt_phase', 0, @isnumeric);
addParameter(func_inputs,'predef_scale_ctrs', [], @isnumeric);
addParameter(func_inputs,'predef_rate_ctrs', [], @isnumeric);
addParameter(func_inputs,'time_const', 1, @isnumeric);
addParameter(func_inputs,'scale_res', 1, @isnumeric);
addParameter(func_inputs,'scale_max', [], @isnumeric);
addParameter(func_inputs,'scale_min', [], @isnumeric);
addParameter(func_inputs,'rate_res', 1, @isnumeric);
addParameter(func_inputs,'rate_max', [], @isnumeric);
addParameter(func_inputs,'rate_min', [], @isnumeric);
parse(func_inputs,varargin{:})

del_cqt_phase = func_inputs.Results.del_cqt_phase;
predef_scale_ctrs = func_inputs.Results.predef_scale_ctrs;
predef_rate_ctrs = func_inputs.Results.predef_rate_ctrs;
time_const = func_inputs.Results.time_const;
scale_res = func_inputs.Results.scale_res;
scale_max = func_inputs.Results.scale_max;
scale_min = func_inputs.Results.scale_min;
rate_res = func_inputs.Results.rate_res;
rate_max = func_inputs.Results.rate_max;
rate_min = func_inputs.Results.rate_min;

comp_scale_ctrs = 1;
comp_rate_ctrs = 1;

if ~isempty(predef_scale_ctrs), comp_scale_ctrs = 0; end
if ~isempty(predef_rate_ctrs), comp_rate_ctrs = 0; end


%% CQT parameters:

fs = cqt_params_in.fs;   % sample rate of the time signal
fmin = cqt_params_in.fmin; % min frequency of the cqt filter bank 
fmax = cqt_params_in.fmax; % max frequency of the cqt filter bank
fres = cqt_params_in.fres; % number of filters per octave
gamma = cqt_params_in.gamma; % linear-Q factor

%% Time-domain signal to CQT

sig_cq_struct = cqt(signal,fres, fs, fmin, fmax, 'rasterize', 'full',...
    'gamma',gamma);
sig_cqt = sig_cq_struct.c;
[n_freq,n_time] = size(sig_cqt); % number of frequency channels and time frames

% observation: this part is invertible for just one
% signal (without any manipulation) but fails in separation
% reason: phase is not linear over harmonics in cqt
% it's actually periodic (over time) with different period
% for each overtone. When you unwrap it over time you can
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

cqt_params_out = sig_cq_struct;
cqt_params_out = rmfield(cqt_params_out,'c');
cqt_params_out.n_freq = n_freq;
cqt_params_out.n_time = n_time;

 
%% Parameters of spectro-temporal filters:

sig_dur = length(signal)/fs; % duration of the signal (in sec)

nfft_scale = n_freq; % min number of fft points along the scale axis 
nfft_rate = n_time; % min number of fft points along the rate axis 

spec_samprate = fres; % sampling rate of the spectral filter (in samples per octave)
temp_samprate = floor(n_time/sig_dur); % sampling rate of the temporal filter (in frames per sec)

% set filter scales and rates if not provided
scale_params = struct();
if comp_scale_ctrs
    scale_params.filt_type = 'scale';
    scale_params.filt_res = scale_res;
    scale_params.filt_nfft = nfft_scale;
    scale_params.samprate = spec_samprate;    
    if ~isempty(scale_max), scale_params.ctr_max = scale_max; end
    if ~isempty(scale_min), scale_params.ctr_min = scale_min; end
    
    scale_ctrs = filt_default_centers_refactored(scale_params); 
end

rate_params = struct();
if comp_rate_ctrs
    rate_params.filt_type = 'rate';
    rate_params.filt_res = rate_res;
    rate_params.filt_nfft = nfft_rate;
    rate_params.samprate = temp_samprate;
    if ~isempty(rate_max), rate_params.ctr_max = rate_max; end
    if ~isempty(rate_min), rate_params.ctr_min = rate_min; end
    
    rate_ctrs = filt_default_centers_refactored(rate_params);
end

%% Spectro-temporal filter bank

disp('Computing the filterbank...');

% scale filter parameters
scale_filt_params = struct('scale_ctrs',scale_ctrs,'nfft_scale',nfft_scale,...
    'spec_samprate',spec_samprate);

% rate filter parameters
rate_filt_params = struct('rate_ctrs',rate_ctrs,'nfft_rate',nfft_rate,...
    'temp_samprate',temp_samprate,'time_const',time_const);

%%% the following line is interesting and worth exploring
% aa = abs(X).*exp(1j*abs(X).*angle(X));
if del_cqt_phase
   [~,fbank_sr_domain] = gen_fbank_scale_rate(scale_filt_params,rate_filt_params);
else
   [~,fbank_sr_domain] = gen_fbank_scale_rate(scale_filt_params,rate_filt_params,...
       sig_cqt);
end

%% CQT to MCFT
disp('Computing the transform...');

if del_cqt_phase
   sig_cqt = abs(sig_cqt);
end

mcft_out=cqt_to_mcft_refactored(sig_cqt,fbank_sr_domain);

end



