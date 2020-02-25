function est_sig_cqt = mcft_to_cqt_light(mcft_in,fbank,filt_ctr_posit,...
    filt_range,sig_cqt_size)

% This function reconstructs the time-frequency representation (CQT) of 
% an audio signal through inverse filtering given the subsampled version
% of the MCFT and the scale-rate domain critically sampled filterbank.
% 
% Inputs: 
% mcft_in: cell array containing mcft coefficients (critically sampled or
%          zero-padded)
% fbank: cell array containing the critically sampled scale-rate 
%          filter bank
% filt_ctr_posit: n_scale * n_rate * 2 matrix containig filter centers (low,band,high)
%          in sample # in the scale-rate domain (1st element: scale, 2nd: rate) 
% filt_range: n_scale * n_rate cell array containing the support of the
%          filter along scale and rate axes in the original sampling
%          grid (transform domain sampled at the actual rate)
% sig_cqt_size: 1*2 matrix containing the number of freq. bins and time frames
%
% Output:
% est_cqt: 2d matrix containing the reconstructed CQT
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Dimensions

[n_scale_ctrs, n_rate_ctrs] = size(fbank);

mcft_size = cellfun(@size, mcft_in, 'UniformOutput', 0);
filt_org_size = cellfun(@size, fbank, 'UniformOutput', 0);


%% CQT reconstruction

% Initialize the estimated CQT
sig_sr_sum = zeros(sig_cqt_size);
fbank_sr_sum = zeros(sig_cqt_size);

for i = 1:n_scale_ctrs
    for j = 1:n_rate_ctrs
        
        % Lengths of critically sampled filters
        scale_len_temp = filt_org_size{i,j}(1);
        rate_len_temp = filt_org_size{i,j}(2);
        
        % Size of the mcft 
        mcft_scale_len = mcft_size{i,j}(1);
        mcft_rate_len = mcft_size{i,j}(2);
         
        % Support of filters on the original grid (original sample rate)
        scale_filt_range = filt_range{i,j}{1};
        rate_filt_range = filt_range{i,j}{2};
        
        % Filter center position (in sample #)
        scale_posit_temp = filt_ctr_posit{i,j}(1);
        rate_posit_temp = filt_ctr_posit{i,j}(2);
        
        % Frequency mapping
        scale_displace = scale_posit_temp - ...
            floor(scale_posit_temp/mcft_scale_len) * mcft_scale_len;        
        rate_displace = rate_posit_temp - ...
            floor(rate_posit_temp/mcft_rate_len) * mcft_rate_len;
        
        % Currect mcft in the scale-rate domain
        sig_sr_temp = fft2(mcft_in{i,j});
        
        % Reverse the frequency shifting
        sig_sr_temp = circshift(sig_sr_temp,-scale_displace,1);
        sig_sr_temp = circshift(sig_sr_temp,-rate_displace,2);
        
        % Remove zero-padding (if exists) and flip to [-pi,pi)
        mcft_scale_idx = mode([mcft_scale_len-floor(scale_len_temp/2)+1:mcft_scale_len,...
            1:ceil(scale_len_temp/2)]-1, mcft_scale_len)+1;
        mcft_rate_idx = mode([mcft_rate_len-floor(rate_len_temp/2)+1:mcft_rate_len,...
            1:ceil(rate_len_temp/2)]-1, mcft_rate_len)+1;
        
        sig_sr_temp = sig_sr_temp(mcft_scale_idx, mcft_rate_idx);
        
        % Prepare the original filter
        filt_scale_idx = [scale_len_temp-floor(scale_len_temp/2)+1:scale_len_temp,...
            1:ceil(scale_len_temp/2)];
        filt_rate_idx = [rate_len_temp-floor(rate_len_temp/2)+1:rate_len_temp,...
            1:ceil(rate_len_temp/2)];
        
        sr_filt_temp = fbank{i,j};
        sr_filt_temp = sr_filt_temp(filt_scale_idx, filt_rate_idx);
        
        % Inverse filter, shift, and add 
        filtered_sig_sr_domain = sig_sr_temp .* conj(sr_filt_temp);
        
        sig_sr_sum(scale_filt_range, rate_filt_range) = ...
            sig_sr_sum(scale_filt_range, rate_filt_range) + filtered_sig_sr_domain;
        
        
        % compute the normalization factor  
        fbank_sr_sum(scale_filt_range, rate_filt_range) = ...
            fbank_sr_sum(scale_filt_range, rate_filt_range) + ...
            sr_filt_temp .* conj(sr_filt_temp);
                       
        
    end
end

% Compute the ratio sig_sr_sum/fbank_sr_sum
sig_sr_ratio = sig_sr_sum ./ (fbank_sr_sum + eps);

% Compute est_sig_cqt
est_sig_cqt = ifft2(sig_sr_ratio);

end


