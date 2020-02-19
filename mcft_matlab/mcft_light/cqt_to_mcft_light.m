function  [mcft_out,fbank_shift,fbank_displace] = cqt_to_mcft_light(sig_cqt,...
    fbank_org,filt_range,filt_ctr_posit,filt_zpad_size) 

% This function computes the coefficients of the subsampled version of 
% the MCFT. It receives the frequency-time representation (CQT) of 
% an audio signal (complex in general) and generates the 4-dimensional 
% representation (scale,rate, frequency, time) by 2d filtering based 
% on the cortical part of Chi's auditory model. 
%
% Inputs:
% sig_cqt: 2d matrix contatining the (complex) time_frequency representation 
%       of an audio signal (log scale frequency, e.g. CQT)
% fbank_sr_domain: 4d matrix containing a bank of filters in the scale-rate 
%       domain
%
% Ouput: 
% mcft_out: ??? will probably be a cell array
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

% dimensions
[scale_nfft,rate_nfft] = size(sig_cqt);
[n_scale_ctrs, n_rate_ctrs] = size(fbank_org);

filt_org_size = cellfun(@size, fbank_org, 'UniformOutput', 0);

% 2D-Fourier transform of the time-frequency representation
sig_sr_domain = fft2(sig_cqt,scale_nfft,rate_nfft);

mcft_out = cell(n_scale_ctrs, n_rate_ctrs);
fbank_shift = cell(n_scale_ctrs, n_rate_ctrs);
fbank_displace = cell(n_scale_ctrs, n_rate_ctrs);
  
for i = 1:n_scale_ctrs
    for j = 1:n_rate_ctrs
        
        scale_len_temp = filt_org_size{i,j}(1);
        rate_len_temp = filt_org_size{i,j}(2);
        
        scale_zpad_len = filt_zpad_size{i,j}(1);
        rate_zpad_len = filt_zpad_size{i,j}(2);
        
        scale_posit_temp = filt_ctr_posit{i,j}(1);
        rate_posit_temp = filt_ctr_posit{i,j}(2);
        
        scale_filt_range = filt_range{i,j}{1};
        rate_filt_range = filt_range{i,j}{2};
        
        sr_filt_temp = fbank_org{i,j};
        
        scale_idx = [ceil(scale_len_temp/2)+1:scale_len_temp, ...
            1:ceil(scale_len_temp/2)];
        rate_idx = [ceil(rate_len_temp/2)+1:rate_len_temp, ...
            1:ceil(rate_len_temp/2)];
        
        mcft_temp = zeros(filt_zpad_size{i,j});
        mcft_temp([end-floor(scale_len_temp/2)+1:end, 1:ceil(scale_len_temp/2)],...
            [end-floor(rate_len_temp/2)+1:end, 1:ceil(rate_len_temp/2)]) = ...
            sig_sr_domain(scale_filt_range, rate_filt_range) .* ...
            sr_filt_temp(scale_idx,rate_idx); 
        % Note: this is the original filter not the zero-padded version
                
        scale_displace = scale_posit_temp - ...
            floor(scale_posit_temp/scale_zpad_len) * scale_zpad_len;        
        rate_displace = rate_posit_temp - ...
            floor(rate_posit_temp/rate_zpad_len) * rate_zpad_len;
        
        mcft_temp = circshift(mcft_temp, scale_displace, 1);
        mcft_temp = circshift(mcft_temp, rate_displace, 2);
        
        mcft_out{i,j} = mcft_temp;
        mcft_out{i,j} = ifft2(mcft_temp);

        
        
        % Now generate zero-padded version and shift (for verification)
        sr_filt_zpad_temp = zeros(filt_zpad_size{i,j});
        sr_filt_zpad_temp([end-floor(scale_len_temp/2)+1:end, 1:ceil(scale_len_temp/2)],...
            [end-floor(rate_len_temp/2)+1:end, 1:ceil(rate_len_temp/2)]) = ...
            sr_filt_temp(scale_idx,rate_idx); 
        
        sr_filt_zpad_temp = circshift(sr_filt_zpad_temp, scale_displace, 1);
        sr_filt_zpad_temp = circshift(sr_filt_zpad_temp, rate_displace, 2);
        
        fbank_shift{i,j} = sr_filt_zpad_temp;
        
        fbank_displace{i,j} = [scale_displace,rate_displace];

        
         
    end
    
end


















