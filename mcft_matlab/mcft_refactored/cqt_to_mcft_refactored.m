function  mcft_out = cqt_to_mcft_refactored(sig_cqt,fbank_sr_domain) 

% This function receives the frequency-time representation (CQT) of 
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
% mcft_out: 4d matrix containing MCFT coefficients
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Dimensions

[n_freq,n_time] = size(sig_cqt);
[num_scale_ctrs,num_rate_ctrs,nfft_scale,nfft_rate] = size(fbank_sr_domain);

%% Time-frequency representation to MCFT:

% 2D-Fourier transform of the time-frequency representation
sig_sr_domain = fft2(sig_cqt,nfft_scale,nfft_rate);

mcft_out = zeros(num_scale_ctrs, num_rate_ctrs, nfft_scale, nfft_rate);

for i = 1:num_scale_ctrs
   for j = 1:num_rate_ctrs
      filt_sr_temp = squeeze(fbank_sr_domain(i,j,:,:));
      
      % filter the signal in scale-rate domain
      filtered_sig_sr_domain = sig_sr_domain .* filt_sr_temp; 
      
      % convert back to the frequency-time domain
      filtered_sig_ft_domain = ifft2(filtered_sig_sr_domain); 
      
      mcft_out(i,j,:,:) = filtered_sig_ft_domain; %(1:n_freq,1:n_time); % remove the zero padding
           
   end
end

end


