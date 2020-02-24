function [fbank,filt_ctr_posit,filt_range] = gen_2d_fbank_light(fbank_params,varargin)

% This function generates a multirate scale-rate (2D) filterbank.
% Bandpass filters are critically sampled by default. If specified, the
% filters can be slightly zeropadded to fit in a matrix format.
%
% Inputs:
% fbank_params: structure array containing scale & rate filter parameters
%       ctr_min: 1*2 vector containing the center of the lowest bandpass 
%                filter in the scale-rate domain [scale_ctr_min,rate_ctr_min]
%       ctr_max: 1*2 vector containing the center of the highest bandpass
%                filter in the scale-rate domain [scale_ctr_max,rate_ctr_max]
%       filt_res: 1*2 vector containing the number of filters (bins) per
%                octave [scale_filt_res,rate_filt_res]
%       nfft: 1*2 vector containing the number of fft points along scale
%             and rate axes [scale_nfft,rate_nfft]
%       samprate: 1*2 vector containing the sample rate along the frequency
%                 axis (cyc/oct) and time axis (cyc/sec)
%                 [spec_samprate,temp_samprate]
%
% Optional input arguments can be provided like this:
% 
%    gen_2d_fbank_light(fbank_params,'min_filt_len',min_filt_len)
%
% The optional arguments must be names(character strings) followed by values:
%
% 'zpadding': 'on' or 'off' (off means critically sampled), default = 'off'
% 'fbank_format': character string, specifies the data structure in which
%                  the output filterbank is stored. The format can be:
%                  - 'cell': for critically sampled and zeropadded cases
%                  - 'matrix': concatenates all cell elements, for critically
%                            sampled and zeropadded cases
%                  - 'struct': containing 4d matrices of same-size bandpass 
%                              or highpass sections of the transform, only 
%                              for the zeropadded case
% 
% 'scale_filt_name': character string specifying the name of the scale filter 
%                    function, default: 'gabor_fourier'
% 'rate_filt_name': character string specifying the name of the scale filter 
%                   function, default: 'gammatone_fourier'
%
% 'min_filt_len': 1*2 vector containing minimum scale and rate filter lengths, 
%                 default: [4,4] (samples)
% 'bw_offset': 1*2 vector containing bandwidth offsets (scale and rate)
%              If bw_offset = 0 the filterbank is constant-Q.
%              bandwidth = 1/Q_fator * ctr + bw_offset 
%              (Q is determined by #bins per octave)
%              bw_offset > 0 improves the resolution of lower filters 
%
% 'sigma': Gabor filter parameter, similar to standard diviation 
% 'time_const': time constant of the exponential term of the gammatone
%               filter
% 
% Outputs: 
% fbank: n_scale * n_rate cell array containing constant-Q/variabale-Q 
%        scale-rate filters
% ctr_posit: n_scale * n_rate * 2 matrix containig filter centers (low,band,high)
%            in sample # in the scale-rate domain (1st element: scale, 2nd: rate) 
% filt_range: n_scale * n_rate cell array containing the support of the
%             filter along scale and rate axes in the original sampling
%             grid (transform domain sampled at the actual rate)
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

% extract filter parameters
ctr_min = fbank_params.ctr_min;
ctr_max = fbank_params.ctr_max;
filt_res = fbank_params.filt_res;
nfft = fbank_params.nfft;
samprate = fbank_params.samprate;

% extract optional parameter values
func_inputs = inputParser;
addParameter(func_inputs,'zpadding', 'off', @ischar);
% addParameter(func_inputs,'fbank_format','cell', @isnumeric);
addParameter(func_inputs,'scale_filt_name','gabor_fourier',@ischar);
addParameter(func_inputs,'rate_filt_name','gammatone_fourier',@ischar);
addParameter(func_inputs,'min_filt_len',[4,4],@isnumeric);
addParameter(func_inputs,'bw_offset',[0,0],@isnumeric);
addParameter(func_inputs,'sigma',1,@isnumeric);
addParameter(func_inputs,'time_const',1,@isnumeric);

parse(func_inputs,varargin{:})
zpadding = func_inputs.Results.zpadding;
% fbank_format = func_inputs.Results.fbank_format;
scale_filt_name = func_inputs.Results.scale_filt_name;
rate_filt_name = func_inputs.Results.rate_filt_name;
min_filt_len = func_inputs.Results.min_filt_len;
bw_offset = func_inputs.Results.bw_offset;
sigma = func_inputs.Results.sigma;
time_const = func_inputs.Results.time_const;

% compute scale filters
[scale_fbank,scale_ctrs,scale_ctr_posit,scale_filt_len] = ...
    gen_1d_fbank_light(scale_filt_name,ctr_min(1),ctr_max(1),filt_res(1),...
    samprate(1),nfft(1),'min_filt_len',min_filt_len(1),'bw_offset',...
    bw_offset(1),'sigma',sigma);

% compute rate filters
[rate_fbank,rate_ctrs,rate_ctr_posit,rate_filt_len] = ...
    gen_1d_fbank_light(rate_filt_name,ctr_min(2),ctr_max(2),filt_res(2),...
    samprate(2),nfft(2),'min_filt_len',min_filt_len(2),'bw_offset',...
    bw_offset(2),'time_const',time_const);

% total number of filters
n_scale_ctrs = length(scale_ctrs);
n_rate_ctrs = length(rate_ctrs);

% add one more high pass filter so that highpass filters can be split into
% up-ward and down-ward later on 
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

% combine scale_filt and rate_filt lengths into one matrix

% 2D filter sizes with and without zero padding
filt_size = zeros(n_scale_ctrs + 1, n_rate_ctrs + 1, 2);

if strcmp(zpadding,'off')    
    filt_size(:,:,1) = repmat(scale_filt_len,1,n_rate_ctrs+1);
    filt_size(:,:,2) = repmat(rate_filt_len.',n_scale_ctrs+1,1);
    
else
    % maximum length bandpass filter
    [scale_bpass_max_len, ~] = max(scale_filt_len(1:n_scale_ctrs/2));
    [rate_bpass_max_len, ~] = max(rate_filt_len(1:n_rate_ctrs/2));

    % highpass lengths
    scale_highpass_len = scale_filt_len(n_scale_ctrs/2+1);
    rate_highpass_len = rate_filt_len(n_rate_ctrs/2+1);
    
    filt_size(:,:,1) = scale_bpass_max_len;
    filt_size(n_scale_ctrs/2 + 1,:,1) = scale_highpass_len; % upward
    filt_size(n_scale_ctrs/2 + 2,:,1) = scale_highpass_len; % downward

    filt_size(:,:,2) = rate_bpass_max_len;
    filt_size(:,n_rate_ctrs/2 + 1,2) = rate_highpass_len; % upward
    filt_size(:,n_rate_ctrs/2 + 2,2) = rate_highpass_len; % downward
    
end

% generate 2d filters
% highpass filters will be split into up-/down-ward filter

fbank = cell(n_scale_ctrs + 1, n_rate_ctrs + 1);

filt_ctr_posit = cell(n_scale_ctrs + 1, n_rate_ctrs + 1);
filt_range = cell(n_scale_ctrs + 1, n_rate_ctrs + 1);

for i = 1:n_scale_ctrs + 1
    for j = 1:n_rate_ctrs + 1
        
    slen_org_temp = scale_filt_len(i);
    sposit_temp = scale_ctr_posit(i);
            
    scale_idx = [ceil(slen_org_temp/2)+1:slen_org_temp, 1:ceil(slen_org_temp/2)];
    scale_filt_range = mod(sposit_temp+(-floor(slen_org_temp/2):ceil(slen_org_temp/2)-1),...
        nfft(1))+1;
       
    rlen_org_temp = rate_filt_len(j);
    rposit_temp = rate_ctr_posit(j);       
        
    rate_idx = [ceil(rlen_org_temp/2)+1:rlen_org_temp, 1:ceil(rlen_org_temp/2)];
    rate_filt_range = mod(rposit_temp+(-floor(rlen_org_temp/2):ceil(rlen_org_temp/2)-1),...
            nfft(2))+1;
                      
    scale_new_len = filt_size(i,j,1);
    sfilt_temp = zeros(scale_new_len,1);
    sfilt_temp([end-floor(slen_org_temp/2)+1:end, 1:ceil(slen_org_temp/2)],:) = ...
                   scale_fbank{i}(scale_idx);
               
    rate_new_len = filt_size(i,j,2);
    rfilt_temp = zeros(rate_new_len,1);
    rfilt_temp([end-floor(rlen_org_temp/2)+1:end, 1:ceil(rlen_org_temp/2)],:) = ...
                   rate_fbank{j}(rate_idx);
               
    sr_filt_temp = sfilt_temp * rfilt_temp.';
    
    fbank{i,j} = sr_filt_temp;
    filt_range{i,j} = {scale_filt_range, rate_filt_range};
    filt_ctr_posit{i,j} = [sposit_temp, rposit_temp];
       
    % separate up-/dow-ward highpass components
    slen_high_half = floor(scale_new_len/2);
    rlen_high_half = floor(rate_new_len/2);
        
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
        sr_filt_half1(2:slen_high_half+1, 2:rlen_high_half+1) = 0;
        sr_filt_half1(slen_high_half+2:end, rlen_high_half+2:end) = 0; 
        
        fbank{i,j} = sr_filt_half1;
        
    end
        
    if condition2
        
        sr_filt_half2 = sr_filt_temp;
        sr_filt_half2(2:slen_high_half+1, rlen_high_half+2:end) = 0;
        sr_filt_half2(slen_high_half+2:end, 2:rlen_high_half+1) = 0;
        
        fbank{i,j} = sr_filt_half2;
         
    end
    

    end
end



end

