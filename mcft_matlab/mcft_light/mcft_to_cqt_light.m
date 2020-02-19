function est_sig_cqt = mcft_to_cqt_light(mcft_in,fbank_org,filt_range,...
    fbank_displace,sig_cqt_size)

% This function reconstructs the time-frequency representation (CQT) of 
% an audio signal through inverse filtering given the subsampled version
% of the MCFT and the scale-rate domain filterbank.
% 
% Inputs: 
% mcft_in: 4d matrix containing the MCFT coefficients
% fbank_sr_domain: 4d matrix containing a bank of scale-rate filters
%
% Output:
% est_cqt: 2d matrix containing the reconstructed time-frequency
%        representation 
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

% dimensions
[n_scale_ctrs, n_rate_ctrs] = size(fbank_org);

filt_zpad_size = cellfun(@size, mcft_in, 'UniformOutput', 0);
filt_org_size = cellfun(@size, fbank_org, 'UniformOutput', 0);


% initialize the estimated cqt
sig_sr_sum = zeros(sig_cqt_size);
fbank_sr_sum = zeros(sig_cqt_size);

for i = 1:n_scale_ctrs
    for j = 1:n_rate_ctrs
        
        scale_len_temp = filt_org_size{i,j}(1);
        rate_len_temp = filt_org_size{i,j}(2);
        
        scale_zpad_len = filt_zpad_size{i,j}(1);
        rate_zpad_len = filt_zpad_size{i,j}(2);
                
        scale_filt_range = filt_range{i,j}{1};
        rate_filt_range = filt_range{i,j}{2};
        
        scale_displace_temp = fbank_displace{i,j}(1);
        rate_displace_temp = fbank_displace{i,j}(2);
       
        sig_sr_temp = fft2(mcft_in{i,j});
        
        % reverse the frequency shifting
        sig_sr_temp = circshift(sig_sr_temp,-scale_displace_temp,1);
        sig_sr_temp = circshift(sig_sr_temp,-rate_displace_temp,2);
        
        % remove zero-padding and flip to [-pi,pi)
        mcft_scale_idx = mode([scale_zpad_len-floor(scale_len_temp/2)+1:scale_zpad_len,...
            1:ceil(scale_len_temp/2)]-1, scale_zpad_len)+1;
        mcft_rate_idx = mode([rate_zpad_len-floor(rate_len_temp/2)+1:rate_zpad_len,...
            1:ceil(rate_len_temp/2)]-1, rate_zpad_len)+1;
        
        sig_sr_temp = sig_sr_temp(mcft_scale_idx, mcft_rate_idx);
        
        
        % prepare the original filter
        filt_scale_idx = [scale_len_temp-floor(scale_len_temp/2)+1:scale_len_temp,...
            1:ceil(scale_len_temp/2)];
        filt_rate_idx = [rate_len_temp-floor(rate_len_temp/2)+1:rate_len_temp,...
            1:ceil(rate_len_temp/2)];
        
        sr_filt_temp = fbank_org{i,j};
        sr_filt_temp = sr_filt_temp(filt_scale_idx, filt_rate_idx);
        
        % overlap and add
        filtered_sig_sr_domain = sig_sr_temp .* conj(sr_filt_temp);
        
        sig_sr_sum(scale_filt_range, rate_filt_range) = ...
            sig_sr_sum(scale_filt_range, rate_filt_range) + filtered_sig_sr_domain;
        
        
        % compute the normalization factor  
        fbank_sr_sum(scale_filt_range, rate_filt_range) = ...
            fbank_sr_sum(scale_filt_range, rate_filt_range) + ...
            sr_filt_temp .* conj(sr_filt_temp);
                       
        
    end
end

% compute the ratio sig_sr_sum/fbank_sr_sum
sig_sr_ratio = sig_sr_sum ./ (fbank_sr_sum + eps);

% compute est_sig_cqt
est_sig_cqt = ifft2(sig_sr_ratio);

end

