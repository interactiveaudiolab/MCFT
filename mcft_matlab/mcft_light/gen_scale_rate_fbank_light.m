function [fbank_org,fbank_zpad,filt_range,filt_ctr_posit,filt_ctrs] = ...
    gen_scale_rate_fbank_light(scale_filt_params,rate_filt_params,varargin)

% This function generates a multirate scale-rate (2D) filterbank.
% Bandpass filters are critically sampled.
%
% Inputs:
% scale_filt_params: structure array containing the parameters of scale filters 
%       ctr_min: center of the lowest bandpass scale filter
%       ctr_max: center of the highest bandpass scale filter
%       filt_res: number of bins (filters) per octave
%       nfft: number of scale fft points
%       spec_samprate: sample rate along the freq. axis (cyc/oct)
%
% rate_filt_params: structure array containing the parameters of rate filters
%       ctr_min: center of the lowest bandpass rate filter
%       ctr_max: center of the highest bandpass rate filter
%       filt_res: number of bins (filters) per octave
%       nfft: number of rate fft points
%       temp_samprate: sample rate along the time axis
%
% Optional input arguments can be provided like this:
% 
%    gen_scale_rate_fbank_light(scale_filt_params,rate_filt_params,...
%                      'min_filt_len',min_filt_len)
%
% The optional arguments much be names(character strings) followed by values:
% 
% 'scale_filt_name': character string specifying the name of the scale filter function
%                    default: 'gabor_fourier'
% 'scale_min_filt_len': minimum scale filter length 
%                    default: 4
% 'scale_bw_offset': bandwidth offset. If bw_offset = 0 the filterbank is constant-Q.
%        bandwidth = 1/Q_fator * ctr + bw_offset (Q is determined by #bins per
%        octave)
%        bw_offset > 0 improves the resolution of lower filters 
% 'sigma': Gabor filter parameter, similar to standard diviation 
%
% 'rate_filt_name': character string specifying the name of the scale filter function
%                    default: 'gammatone_fourier'
% 'rate_min_filt_len': minimum scale filter length 
%                    default: 4
% 'rate_bw_offset': bandwidth offset. If bw_offset = 0 the filterbank is constant-Q.
%        bandwidth = 1/Q_fator * ctr + bw_offset (Q is determined by #bins per
%        octave)
%        bw_offset > 0 improves the resolution of lower filters 
% 'time_const': time constant of the exponential term
% 
% Outputs: 
% fbank: n_scale * n_rate cell array of containing constant-Q/variabale-Q 
%        scale-rate filters
% ctr_shift: n_scale * n_rate * 2 matrix containig shifts between the center 
%            frequencies in the scale-rate domain (1st element: scale, 2nd: rate)            
% filt_size: n_scale * n_rate * 2 matrix containing filter sizes
%

% Extract scale parameters
scale_ctr_min = scale_filt_params.ctr_min;
scale_ctr_max = scale_filt_params.ctr_max;
scale_filt_res = scale_filt_params.filt_res;
scale_nfft = scale_filt_params.nfft;
spec_samprate = scale_filt_params.samprate;

% Extract rate parameters
rate_ctr_min = rate_filt_params.ctr_min;
rate_ctr_max = rate_filt_params.ctr_max;
rate_filt_res = rate_filt_params.filt_res;
rate_nfft = rate_filt_params.nfft;
temp_samprate = rate_filt_params.samprate;

% Extract optional parameter values
func_inputs = inputParser;
addParameter(func_inputs,'scale_filt_name','gabor_fourier',@ischar);
addParameter(func_inputs,'scale_min_filt_len',4,@isnumeric);
addParameter(func_inputs,'scale_bw_offset',0,@isnumeric);
addParameter(func_inputs,'sigma',1,@isnumeric);

addParameter(func_inputs,'rate_filt_name','gammatone_fourier',@ischar);
addParameter(func_inputs,'rate_min_filt_len',4,@isnumeric);
addParameter(func_inputs,'rate_bw_offset',0,@isnumeric);
addParameter(func_inputs,'time_const',1,@isnumeric);

parse(func_inputs,varargin{:})
scale_filt_name = func_inputs.Results.scale_filt_name;
scale_min_filt_len = func_inputs.Results.scale_min_filt_len;
scale_bw_offset = func_inputs.Results.scale_bw_offset;
sigma = func_inputs.Results.sigma;

rate_filt_name = func_inputs.Results.rate_filt_name;
rate_min_filt_len = func_inputs.Results.rate_min_filt_len;
rate_bw_offset = func_inputs.Results.rate_bw_offset;
time_const = func_inputs.Results.time_const;

% Compute scale filters
[scale_fbank,scale_ctrs,scale_ctr_posit,scale_filt_len] = ...
    gen_scale_fbank_light(scale_ctr_min,scale_ctr_max,scale_filt_res,...
    spec_samprate,scale_nfft,'filt_name',scale_filt_name,...
    'min_filt_len',scale_min_filt_len,'bw_offset',scale_bw_offset,'sigma',sigma);
 
% Compute rate filters
[rate_fbank,rate_ctrs,rate_ctr_posit,rate_filt_len] = ...
    gen_rate_fbank_light(rate_ctr_min,rate_ctr_max,rate_filt_res,...
    temp_samprate,rate_nfft,'filt_name',rate_filt_name,...
    'min_filt_len',rate_min_filt_len,'bw_offset',rate_bw_offset,'time_const',time_const);

filt_ctrs = {scale_ctrs, rate_ctrs};

n_scale_ctrs = length(scale_ctrs);
n_rate_ctrs = length(rate_ctrs);

% Maximum length bandpass filter
[scale_bpass_max_len, ~] = max(scale_filt_len(1:n_scale_ctrs/2));
[rate_bpass_max_len, ~] = max(rate_filt_len(1:n_rate_ctrs/2));

% Highpass lengths
scale_highpass_len = scale_filt_len(n_scale_ctrs/2+1);
rate_highpass_len = rate_filt_len(n_rate_ctrs/2+1);

%%%%%%%%% combine scale_filt and rate_filt lengths into one matrix

% 2D filter sizes with zero padding
filt_zpad_size = zeros(n_scale_ctrs + 1, n_rate_ctrs + 1, 2);

filt_zpad_size(:,:,1) = scale_bpass_max_len;
filt_zpad_size(n_scale_ctrs/2 + 1,:,1) = scale_highpass_len; % upward
filt_zpad_size(n_scale_ctrs/2 + 2,:,1) = scale_highpass_len; % downward

filt_zpad_size(:,:,2) = rate_bpass_max_len;
filt_zpad_size(:,n_rate_ctrs/2 + 1,2) = rate_highpass_len; % upward
filt_zpad_size(:,n_rate_ctrs/2 + 2,2) = rate_highpass_len; % downward

%%%%%%%%% generate 2d filters
% highpass filters will be split into up-/down-ward filter
% therefore -> one more set of filters added
scale_fbank = [scale_fbank(1:n_scale_ctrs/2+1);...
               scale_fbank(n_scale_ctrs/2+1);...
               scale_fbank(n_scale_ctrs/2+2:end)];
           
rate_fbank = [rate_fbank(1:n_rate_ctrs/2+1);...
              rate_fbank(n_rate_ctrs/2+1);...
              rate_fbank(n_rate_ctrs/2+2:end)];

scale_filt_len = [scale_filt_len(1:n_scale_ctrs/2+1);...
                  scale_filt_len(n_scale_ctrs/2+1);...
                  scale_filt_len(n_scale_ctrs/2+2:end)];
              
rate_filt_len = [rate_filt_len(1:n_rate_ctrs/2+1);...
                 rate_filt_len(n_rate_ctrs/2+1);...
                 rate_filt_len(n_rate_ctrs/2+2:end)];
             
scale_ctr_posit = [scale_ctr_posit(1:n_scale_ctrs/2+1);...
                   scale_ctr_posit(n_scale_ctrs/2+1);...
                   scale_ctr_posit(n_scale_ctrs/2+2:end)];
               
rate_ctr_posit = [rate_ctr_posit(1:n_rate_ctrs/2+1);...
                  rate_ctr_posit(n_rate_ctrs/2+1);...
                  rate_ctr_posit(n_rate_ctrs/2+2:end)];


fbank_org = cell(n_scale_ctrs + 1, n_rate_ctrs + 1);
fbank_zpad = cell(n_scale_ctrs + 1, n_rate_ctrs + 1);

filt_range = cell(n_scale_ctrs + 1, n_rate_ctrs + 1);
filt_ctr_posit = cell(n_scale_ctrs + 1, n_rate_ctrs + 1);

for i = 1:n_scale_ctrs + 1
    for j = 1:n_rate_ctrs + 1
    
    slen_temp = scale_filt_len(i);
    sposit_temp = scale_ctr_posit(i);
            
    scale_idx = [ceil(slen_temp/2)+1:slen_temp, 1:ceil(slen_temp/2)];
    scale_filt_range = mod(sposit_temp+(-floor(slen_temp/2):ceil(slen_temp/2)-1),...
        scale_nfft)+1;
    
                     
    rlen_temp = rate_filt_len(j);
    rposit_temp = rate_ctr_posit(j);       
        
    rate_idx = [ceil(rlen_temp/2)+1:rlen_temp, 1:ceil(rlen_temp/2)];
    rate_filt_range = mod(rposit_temp+(-floor(rlen_temp/2):ceil(rlen_temp/2)-1),...
            rate_nfft)+1;
        
    sr_filt_temp =  scale_fbank{i} * rate_fbank{j}.';
    fbank_org{i,j} = sr_filt_temp;
              
    scale_zpad_len = filt_zpad_size(i,j,1);
    sfilt_zpad_temp = zeros(scale_zpad_len,1);
    sfilt_zpad_temp([end-floor(slen_temp/2)+1:end, 1:ceil(slen_temp/2)],:) = ...
                   scale_fbank{i}(scale_idx);
               
    rate_zpad_len = filt_zpad_size(i,j,2);
    rfilt_zpad_temp = zeros(rate_zpad_len,1);
    rfilt_zpad_temp([end-floor(rlen_temp/2)+1:end, 1:ceil(rlen_temp/2)],:) = ...
                   rate_fbank{j}(rate_idx);
               
    sr_filt_zpad_temp = sfilt_zpad_temp * rfilt_zpad_temp.';
    
    fbank_zpad{i,j} = sr_filt_zpad_temp;
    filt_range{i,j} = {scale_filt_range, rate_filt_range};
    filt_ctr_posit{i,j} = [sposit_temp, rposit_temp];
       
    % separate up-/dow-ward highpass components
    slen = floor(slen_temp/2);
    rlen = floor(rlen_temp/2);
    
    slen_zpad = floor(scale_zpad_len/2);
    rlen_zpad = floor(rate_zpad_len/2);
    
    condition1 = (i == n_scale_ctrs/2 + 1 && j > 1 && j < n_rate_ctrs/2 + 1) || ...
                 (i == n_scale_ctrs/2 + 2 && j > n_rate_ctrs/2 + 2) || ...
                 (i == n_scale_ctrs/2 + 1 && j ==  n_rate_ctrs/2 + 2) || ...
                 (i == n_scale_ctrs/2 + 2 && j ==  n_rate_ctrs/2 + 1) || ...
                 (j == n_rate_ctrs/2 + 1  && i > 1 && i < n_scale_ctrs/2 + 1) || ...
                 (j == n_rate_ctrs/2 + 2  && i > n_scale_ctrs/2 + 2);
             
    condition2 = (i == n_scale_ctrs/2 + 1 && j > n_rate_ctrs/2 + 2) || ...
                 (i == n_scale_ctrs/2 + 2 && j > 1 && j < n_rate_ctrs/2 + 1) || ...
                 (i == n_scale_ctrs/2 + 1 && j ==  n_rate_ctrs/2 + 1) || ...
                 (i == n_scale_ctrs/2 + 2 && j ==  n_rate_ctrs/2 + 2) || ...
                 (j == n_rate_ctrs/2 + 1  && i > n_scale_ctrs/2 + 2) || ...
                 (j == n_rate_ctrs/2 + 2  && i > 1 && i < n_scale_ctrs/2 + 1);
                 
        
    if condition1
        
        sr_filt_half1 = sr_filt_temp;
        sr_filt_half1(2:slen+1, 2:rlen+1) = 0;
        sr_filt_half1(slen+2:end, rlen+2:end) = 0; 
        
        fbank_org{i,j} = sr_filt_half1;
        
        sr_filt_zpad_half1 = sr_filt_zpad_temp;
        sr_filt_zpad_half1(2:slen_zpad+1, 2:rlen_zpad+1) = 0;
        sr_filt_zpad_half1(slen_zpad+2:end, rlen_zpad+2:end) = 0; 
        
        fbank_zpad{i,j} = sr_filt_zpad_half1;
        
    end
        
    if condition2
        
        sr_filt_half2 = sr_filt_temp;
        sr_filt_half2(2:slen+1, rlen+2:end) = 0;
        sr_filt_half2(slen+2:end, 2:rlen+1) = 0;
        
        fbank_org{i,j} = sr_filt_half2;
        
        sr_filt_zpad_half2 = sr_filt_zpad_temp;
        sr_filt_zpad_half2(2:slen_zpad+1, rlen_zpad+2:end) = 0;
        sr_filt_zpad_half2(slen_zpad+2:end, 2:rlen_zpad+1) = 0; 
        
        fbank_zpad{i,j} = sr_filt_zpad_half2;
         
    end
    

    end
end


end











