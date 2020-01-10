function est_sig_cqt = mcft_to_cqt_refactored(mcft_in,fbank_sr_domain)

% This function reconstructs the time-frequency representation (CQT) of 
% an audio signal through inverse filtering given the 4-dimensioanl
% MCFT representation and the scale-rate domain filterbank.
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


%% Dimensions

[num_scale_ctrs,num_rate_ctrs,n_freq,n_time] = size(mcft_in);
[~,~,nfft_scale,nfft_rate] = size(fbank_sr_domain);

%% MCFT to time-frequency representation

fbank_sr_sum = 0;
sig_sr_sum = 0;
for i = 1:num_scale_ctrs
    for j = 1:num_rate_ctrs
         
         mcft_temp = squeeze(mcft_in(i,j,:,:));
         sig_sr_temp = fft2(mcft_temp,nfft_scale,nfft_rate);     
         filt_sr_temp = squeeze(fbank_sr_domain(i,j,:,:));
         
         filtered_sig_sr_domain = sig_sr_temp .* conj(filt_sr_temp);
                  
         fbank_sr_sum = fbank_sr_sum + filt_sr_temp .* conj(filt_sr_temp);
         sig_sr_sum = sig_sr_sum + filtered_sig_sr_domain;
         
    end
end
                  
% compute the ratio sig_sr_sum/fbank_sr_sum
sig_sr_ratio = sig_sr_sum ./ (fbank_sr_sum + eps);

% compute est_sig_cqt
est_sig_cqt = ifft2(sig_sr_ratio);
est_sig_cqt = est_sig_cqt(1:n_freq, 1:n_time);

end
