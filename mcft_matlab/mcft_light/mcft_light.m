function [mcft_out,inv_bundle] = mcft_light(signal,cqt_params_in,fbank_params,varargin)

% This function computes a subsampled version of the MCFT, where subsampled
% means critically sampled (default) or critically sample + slight
% zeropadding. 
% 
% Inputs:
% signal: vector containing samples of the time-domain signal
%
% cqt_params_in: structure array containing CQT parameters:
%        fs: sampling rate of the audio signal
%        fmin: minimum frequency of analysis
%        fmax: maximum frequency of analysis
%        fres: frequency resolution (# of bins per octave)
%        gamma: linear-Q factor
%
% fbank_params: structure array containing scale & rate filter parameters
%       ctr_min: 1*2 vector containing the center of the lowest bandpass 
%                filter in the scale-rate domain [scale_ctr_min,rate_ctr_min]
%       ctr_max: 1*2 vector containing the center of the highest bandpass
%                filter in the scale-rate domain [scale_ctr_max,rate_ctr_max]
%       filt_res: 1*2 vector containing the number of filters (bins) per
%                octave [scale_filt_res,rate_filt_res]
% 
%
% Optional inputs:
% 'zpadding': 'on' or 'off' (off means critically sampled), default = 'off'
%
% 'output_format': character string, specifies the data structure in which
%                  the mcft coefficients are stored. The format can be:
%                  - 'cell': for critically sampled and zeropadded cases
%                  - 'matrix': concatenates all cell elements, for critically
%                            sampled and zeropadded cases
%                  - 'struct': containing 4d matrices of same-size bandpass 
%                              or highpass sections of the transform, only 
%                              for the zeropadded case
%
% 'del_cqt_phase': boolean, indicating whether the cqt phase 
%                 is to be deleted or included in the 2d filtering process. 
%                 Note: the output representation is invertible (back to 
%                       time domain) only if the phase is preserved (default = 0)
%
% 'sigma': Gabor filter parameter, similar to standard diviation, default = 1
% 'time_const': time constant of the temporal/rate filter, default = 1
% 'min_filt_len': 1*2 vector containing minimum scale and rate filter lengths, 
%                 default: [4,4] (samples)
% 'bw_offset': 1*2 vector containing bandwidth offsets (scale and rate)
%              If bw_offset = 0 the filterbank is constant-Q.
%              bandwidth = 1/Q_fator * ctr + bw_offset 
%              (Q is determined by #bins per octave)
%              bw_offset > 0 improves the resolution of lower filters
%
% 
% Outputs:
% mcft_out: cell array, matrix, or structure array containing mcft
%           coefficients
% inv_bundle: structure array containing the information required 
%             for the inverse transform (to reconstruct the time-domain signal)
% 
%       cqt_params: structure array containing all CQT properties 
%       fbank: cell array containing the critically sampled scale-rate filter bank
%       filt_ctr_posit: n_scale * n_rate * 2 matrix containig filter centers (low,band,high)
%            in sample # in the scale-rate domain (1st element: scale, 2nd: rate) 
%       filt_range: n_scale * n_rate cell array containing the support of the
%             filter along scale and rate axes in the original sampling
%             grid (transform domain sampled at the actual rate)
%       
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Input check

% extract optional inputs
func_inputs = inputParser;
addParameter(func_inputs,'zpadding', 'off', @ischar);
addParameter(func_inputs,'output_format','cell', @ischar);
addParameter(func_inputs,'del_cqt_phase', 0, @isnumeric);
addParameter(func_inputs,'sigma',1, @isnumeric);
addParameter(func_inputs,'time_const',1, @isnumeric);
addParameter(func_inputs,'min_filt_len',[4,4], @isnumeric);
addParameter(func_inputs,'bw_offset',[0,0], @isnumeric);
parse(func_inputs,varargin{:})

zpadding = func_inputs.Results.zpadding;
output_format = func_inputs.Results.output_format;
del_cqt_phase = func_inputs.Results.del_cqt_phase;
sigma = func_inputs.Results.sigma;
time_const = func_inputs.Results.time_const;
min_filt_len = func_inputs.Results.min_filt_len;
bw_offset = func_inputs.Results.bw_offset;


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

cqt_params_out = sig_cq_struct;
cqt_params_out = rmfield(cqt_params_out,'c');
cqt_params_out.n_freq = n_freq;
cqt_params_out.n_time = n_time;

inv_bundle = struct('cqt_params',cqt_params_out);

%% Parameters of spectro-temporal filters:

ctr_min = fbank_params.ctr_min;
ctr_max = fbank_params.ctr_max;
filt_res = fbank_params.filt_res;

sig_dur = length(signal)/fs; % duration of the signal (in sec)

scale_nfft = n_freq; % min number of fft points along the scale axis 
rate_nfft = n_time; % min number of fft points along the rate axis 

spec_samprate = fres; % sampling rate of the spectral filter (in samples per octave)
temp_samprate = floor(n_time/sig_dur); % sampling rate of the temporal filter (in frames per sec)
                      
fbank_params = struct('ctr_min',ctr_min,'ctr_max',ctr_max,'filt_res',filt_res,...
    'nfft',[scale_nfft,rate_nfft],'samprate',[spec_samprate,temp_samprate]);                       
                       
%% Spectro-temporal filter bank

disp('Computing the filterbank...');

% Note: the filters are critically sampled
[fbank,filt_ctr_posit,filt_range] = gen_2d_fbank_light(fbank_params,...
    'min_filt_len',min_filt_len,'bw_offset',bw_offset,'sigma',sigma,...
    'time_const',time_const);

% if del_cqt_phase
%    [~,fbank_sr_domain] = gen_fbank_scale_rate(scale_filt_params,rate_filt_params);
% else
%    [~,fbank_sr_domain] = gen_fbank_scale_rate(scale_filt_params,rate_filt_params,...
%        sig_cqt);
% end

inv_bundle.fbank = fbank;
inv_bundle.filt_ctr_posit = filt_ctr_posit;
inv_bundle.filt_range = filt_range;


%% CQT to MCFT

if del_cqt_phase
   sig_cqt = abs(sig_cqt);
end

disp('Computing the transform...');

mcft_cell = cqt_to_mcft_light(sig_cqt,fbank,filt_range,filt_ctr_posit,...
    'zpadding',zpadding);

%% Output format

switch output_format
    
    case 'cell' % critically sampled or zero-padded
        mcft_out = mcft_cell;
        
    case 'matrix' % critically sampled or zero-padded
        mcft_out = cell2mat(mcft_cell);
        
    case 'struct' % only zero-padded
        mcft_out = cell_to_struct(mcft_cell,fbank);  
end

      
end 

function mcft_struct = cell_to_struct(mcft_cell,fbank)

% This function groups the same-size sections of the filtered CQT and
% stores them in 4D matrices, which are in trun stored in a structure array

% Note: lowpass filters are zero-padded to the bandpass size 

% Number of filter centers
[n_scale_ctrs, n_rate_ctrs] = size(fbank);

% Size of zero-padded slices
mcft_zpad_size = cellfun(@size, mcft_cell, 'UniformOutput',0);

% Highpass filter indices
scale_high_idx = ceil(n_scale_ctrs/2);
rate_high_idx = ceil(n_rate_ctrs/2);

% Same-size groups
sband_rband_size = mcft_zpad_size{2,2}; 
shigh_rband_size = mcft_zpad_size{scale_high_idx,2};
shigh_rhigh_size = mcft_zpad_size{scale_high_idx,rate_high_idx};
sband_rhigh_size = mcft_zpad_size{2,rate_high_idx};

% All bandpass- and lowpass-filtered slices
mcft_band_cell = [mcft_cell(1:scale_high_idx-1, 1:rate_high_idx-1),...
                  mcft_cell(1:scale_high_idx-1, rate_high_idx+2:end);...
                  mcft_cell(scale_high_idx+2:end, 1:rate_high_idx-1),...
                  mcft_cell(scale_high_idx+2:end, rate_high_idx+2:end)];
              
mcft_band = zeros(n_scale_ctrs-2, n_rate_ctrs-2, ...
    sband_rband_size(1),sband_rband_size(2));

for i = 1:n_scale_ctrs-2
    for j = 1:n_rate_ctrs-2
        mcft_band(i,j,:,:) = mcft_band_cell{i,j};
    end
end
     
% All highpass-scale-/bandpass-rate-filterered slices
mcft_shigh_rband_cell = [mcft_cell(scale_high_idx:scale_high_idx+1,1:rate_high_idx-1),...
                    mcft_cell(scale_high_idx:scale_high_idx+1,rate_high_idx+2:end)];
                  
mcft_shigh_rband = zeros(2, n_rate_ctrs-2, shigh_rband_size(1),shigh_rband_size(2));

for j = 1:n_rate_ctrs-2
    mcft_shigh_rband(1,j,:,:) = mcft_shigh_rband_cell{1,j};
    mcft_shigh_rband(2,j,:,:) = mcft_shigh_rband_cell{2,j};
end    

% All bandpass-scale-/highpass-rate-filtered slices
mcft_sband_rhigh_cell = [mcft_cell(1:scale_high_idx-1,rate_high_idx:rate_high_idx+1);...
                    mcft_cell(scale_high_idx+2:end,rate_high_idx:rate_high_idx+1)];
                  
mcft_sband_rhigh = zeros(n_scale_ctrs-2,2,sband_rhigh_size(1),sband_rhigh_size(2));

for i = 1:n_scale_ctrs-2
    mcft_sband_rhigh(i,1,:,:) = mcft_sband_rhigh_cell{i,1};
    mcft_sband_rhigh(i,2,:,:) = mcft_sband_rhigh_cell{i,2};
end    

% All highpass-scale-/highpass-rate-filtered slices
mcft_shigh_rhigh_cell = mcft_cell(scale_high_idx:scale_high_idx+1,rate_high_idx:rate_high_idx+1);

mcft_shigh_rhigh = zeros(2,2,shigh_rhigh_size(1),shigh_rhigh_size(2));
for i = 1:2
    mcft_shigh_rhigh(i,1,:,:) = mcft_shigh_rhigh_cell{i,1};
    mcft_shigh_rhigh(i,2,:,:) = mcft_shigh_rhigh_cell{i,2};
end

% Output structure array
mcft_struct = struct('sband_rband',mcft_band,'shight_rband',mcft_shigh_rband,...
                         'sband_rhigh',mcft_sband_rhigh,'shigh_rhigh',mcft_shigh_rhigh);
                    
end

