function  [mcft_out,fbank_displace] = ...
    cqt_to_mcft_light(sig_cqt,fbank,filt_range,filt_ctr_posit,varargin) 

% This function computes the coefficients of the subsampled version of 
% the MCFT. It receives the frequency-time representation (CQT) of 
% an audio signal (complex in general) and generates the 4-dimensional 
% representation (scale,rate, frequency, time) by 2d filtering based 
% on the cortical part of Chi's auditory model. 
%
% Inputs:
% sig_cqt: 2d matrix contatining the (complex) time_frequency representation 
%       of an audio signal (log scale frequency, e.g. CQT)
% fbank: cell array containing a bank of filters in the scale-rate domain
%        note: all the filters are critically sampled        
% filt_range: n_scale * n_rate cell array containing the support of the
%             filter along scale and rate axes in the original sampling
%             grid (transform domain sampled at the actual rate)
% filt_ctr_posit: n_scale * n_rate * 2 matrix containig filter centers (low,band,high)
%            in sample # in the scale-rate domain (1st element: scale, 2nd: rate)  
%
% Optional input arguments can be provided like this:
% 
%    cqt_to_mcft(sig_cqt,fbank,filt_range,filt_ctr_posit,'zpadding',zpadding)
%
% The optional arguments must be names(character strings) followed by values:
%
% 'zpadding': 'on' or 'off' (off means critically sampled), default = 'off'
% 
% Ouput: 
% mcft_out: n_scale * n_rate cell array containing the MCFT coefficients 
% fbank_displace: n_scale * n_rate cell array constaining the filter shifts
%                 which are used in frequency mapping (for phase correction)
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

% extract optional parameter values
func_inputs = inputParser;
addParameter(func_inputs,'zpadding', 'off', @ischar);

parse(func_inputs,varargin{:})
zpadding = func_inputs.Results.zpadding;

% dimensions
[scale_nfft,rate_nfft] = size(sig_cqt);
[n_scale_ctrs, n_rate_ctrs] = size(fbank);

filt_org_size = cellfun(@size, fbank, 'UniformOutput', 0);

% 2D filter sizes with zero padding
if strcmp(zpadding,'on')    
    filt_zpad_size = zeros(n_scale_ctrs, n_rate_ctrs, 2);
    
    scale_filt_len = cell2mat(filt_org_size(:,1));
    scale_filt_len = scale_filt_len(:,1);
    
    rate_filt_len = cell2mat(filt_org_size(1,:).');
    rate_filt_len = rate_filt_len(:,2);

    % maximum length bandpass filter
    [scale_bpass_max_len, ~] = max(scale_filt_len(1:floor(n_scale_ctrs/2)));
    [rate_bpass_max_len, ~] = max(rate_filt_len(1:floor(n_rate_ctrs/2)));

    % highpass lengths
    scale_highpass_len = scale_filt_len(floor(n_scale_ctrs/2)+1);
    rate_highpass_len = rate_filt_len(floor(n_rate_ctrs/2)+1);
    
    filt_zpad_size(:,:,1) = scale_bpass_max_len;
    filt_zpad_size(floor(n_scale_ctrs/2) + 1,:,1) = scale_highpass_len; % upward
    filt_zpad_size(floor(n_scale_ctrs/2) + 2,:,1) = scale_highpass_len; % downward

    filt_zpad_size(:,:,2) = rate_bpass_max_len;
    filt_zpad_size(:,floor(n_rate_ctrs/2) + 1,2) = rate_highpass_len; % upward
    filt_zpad_size(:,floor(n_rate_ctrs/2) + 2,2) = rate_highpass_len; % downward
    
end

% 2D-Fourier transform of the time-frequency representation
sig_sr_domain = fft2(sig_cqt,scale_nfft,rate_nfft);

mcft_out = cell(n_scale_ctrs, n_rate_ctrs);
fbank_displace = cell(n_scale_ctrs, n_rate_ctrs);
  
for i = 1:n_scale_ctrs
    for j = 1:n_rate_ctrs
        
        % lengths of critically sampled filters
        scale_len_temp = filt_org_size{i,j}(1);
        rate_len_temp = filt_org_size{i,j}(2);
        
        % length of the mcft (critically sampled or zero-padded)
        if strcmp(zpadding,'off')
            scale_new_len = scale_len_temp;
            rate_new_len = rate_len_temp;
        else
            scale_new_len = filt_zpad_size(i,j,1);
            rate_new_len = filt_zpad_size(i,j,2);
        end
        
        % filter center position (in sample #)
        scale_posit_temp = filt_ctr_posit{i,j}(1);
        rate_posit_temp = filt_ctr_posit{i,j}(2);
        
        % support of filters on the original grid (original sample rate)
        scale_filt_range = filt_range{i,j}{1};
        rate_filt_range = filt_range{i,j}{2};
        
        % critically sampled filter
        sr_filt_temp = fbank{i,j};
        
        % indices of the filter samples (-pi,pi]
        scale_idx = [ceil(scale_len_temp/2)+1:scale_len_temp, ...
            1:ceil(scale_len_temp/2)];
        rate_idx = [ceil(rate_len_temp/2)+1:rate_len_temp, ...
            1:ceil(rate_len_temp/2)];
        
        % compute mcft coefficients
        mcft_temp = zeros(scale_new_len,rate_new_len);
        mcft_temp([end-floor(scale_len_temp/2)+1:end, 1:ceil(scale_len_temp/2)],...
            [end-floor(rate_len_temp/2)+1:end, 1:ceil(rate_len_temp/2)]) = ...
            sig_sr_domain(scale_filt_range, rate_filt_range) .* ...
            sr_filt_temp(scale_idx,rate_idx); 
        % Note: this is the original filter not the zero-padded version
                
        scale_displace = scale_posit_temp - ...
            floor(scale_posit_temp/scale_new_len) * scale_new_len;        
        rate_displace = rate_posit_temp - ...
            floor(rate_posit_temp/rate_new_len) * rate_new_len;
        
        mcft_temp = circshift(mcft_temp, scale_displace, 1);
        mcft_temp = circshift(mcft_temp, rate_displace, 2);
        
        mcft_out{i,j} = ifft2(mcft_temp);
        
        fbank_displace{i,j} = [scale_displace,rate_displace];

        
         
    end
end


end














