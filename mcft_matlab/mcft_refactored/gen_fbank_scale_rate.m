function [fbank_tf_domain,fbank_sr_domain] = gen_fbank_scale_rate(scale_filt_params,...
    rate_filt_params, complex_specgram)

% This function generates a scale-rate-domain bank of up-/down-ward,
% filters. The filterbank will be tuned to the passband of a target 
% signal if specified.
% 
% Inputs: 
% scale_filt_params: structure array containing the parameters of scale filters 
%       scale_ctrs: vector containing filter centers along the scale axis
%       nfft_scale: number of scale fft points
%       spec_samprate: sample rate along the freq. axis (cyc/oct)
% 
% rate_filt_params: structure array containing the parameters of rate filters
%       rate_ctrs: vector containing filter centers along the rate axis
%       nfft_rate: number of rate fft points
%       temp_samprate: sample rate along the time axis
%       time_const: exponent coefficient of the exponential term
% 
% complex_specgram (optional): n_freq by n_time matrix containing a 2d 
%       complex spectrogram. 
%       If provided, the function will return a filter bank that is modulated 
%       with the phase of the spectrogram. Otherwise, the function will return 
%       the original set of filters.   
% 
% Outputs: 
% fbank_tf_domain: n_scale * (2*n_rate) * nfft_scale * nfft_rate matix 
%       containing the time-freq-domain filter bank
% fbank_sr_domain: n_scale * (2*n_rate) * nfft_scale * nfft_rate matrix 
%       containing the scale-rate-domain filter bank
% Note: nfft_scale >= n_freq and nfft_rate >= n_time, where 
%       [n_freq , n_time] = size(complex_specgram)
%
% Note: the first and last filters in scale and rate ranges are assumed
%       to be lowpass and highpass respectively
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)


%% Input check
mod_filter = 0;
if nargin == 3
    mod_filter = 1;
end

%% Extract filter parameters

% scale filter parameters
scale_ctrs = scale_filt_params.scale_ctrs;
nfft_scale = scale_filt_params.nfft_scale;

% rate filter parameters
rate_ctrs = rate_filt_params.rate_ctrs;
nfft_rate = rate_filt_params.nfft_rate;

%% Parameters and dimensions

nfft_scale = nfft_scale + mod(nfft_scale,2); % set nfft to the next even number
nfft_rate = nfft_rate + mod(nfft_rate,2);

scale_params = rmfield(scale_filt_params, {'scale_ctrs','nfft_scale'});
scale_params.spec_filt_len = nfft_scale;

rate_params = rmfield(rate_filt_params, {'rate_ctrs','nfft_rate'});
rate_params.temp_filt_len = nfft_rate;

num_scale_ctrs = length(scale_ctrs);
num_rate_ctrs = length(rate_ctrs);

%% Generate the filterbank

% Filter modulation factor (can be thought as a pre-filtering stage)
if mod_filter     
   sig_ft = ifft2(fft2(complex_specgram,nfft_scale,nfft_rate)); % for dimension adjustment
   sig_ft_ph = angle(sig_ft); % complete phase
   % uncomment this line to compare to the python code
   % otherwise fft phase difference results in large error values
   % sig_ft_ph = angle(sig_ft).*(abs(sig_ft) > 1e-10); % includes phase of nonzero elements only
   filt_phase_factor = exp(1j .* sig_ft_ph);
   
else
   filt_phase_factor = 1;
   
end

fbank_tf_domain = zeros(num_scale_ctrs, 2*num_rate_ctrs, nfft_scale, nfft_rate);
fbank_sr_domain = zeros(num_scale_ctrs, 2*num_rate_ctrs, nfft_scale, nfft_rate);
for i = 1:num_scale_ctrs
    if i == 1 %scale_ctrs(i) == scale_ctrs(1) 
        scale_params.spec_filt_type = 'lowpass'; 
    elseif i == num_scale_ctrs %scale_ctrs(i) == scale_ctrs(end)
        scale_params.spec_filt_type = 'highpass'; 
    else
        scale_params.spec_filt_type = 'bandpass';
    end
        
    for j = 1:num_rate_ctrs
         if j == 1 %rate_ctrs(j) == rate_ctrs(1)
           rate_params.temp_filt_type = 'lowpass';
         elseif j == num_rate_ctrs %rate_ctrs(j) == rate_ctrs(end)
           rate_params.temp_filt_type = 'highpass';
         else
           rate_params.temp_filt_type = 'bandpass';
         end
                  
         % generate two analytic filters (one upward and one downward)
         % for the current (Scale,Rate) values
         
         % upward
         [filt_tf_up,~] = gen_filt_scale_rate(scale_ctrs(i),rate_ctrs(j),...
             scale_params,rate_params,'up');         
         filt_tf_up_tuned = filt_tf_up .* filt_phase_factor;
         fbank_tf_domain(i,num_rate_ctrs-j+1,:,:) = filt_tf_up_tuned;         
         filt_sr_up = fft2(filt_tf_up_tuned);
         fbank_sr_domain(i,num_rate_ctrs-j+1,:,:) = filt_sr_up;         
                  
         % downward
         [filt_tf_down,~] = gen_filt_scale_rate(scale_ctrs(i),rate_ctrs(j),...
             scale_params,rate_params,'down');         
         filt_tf_down_tuned = filt_tf_down .* filt_phase_factor;         
         fbank_tf_domain(i,num_rate_ctrs+j,:,:) = filt_tf_down_tuned;
                           
         filt_sr_down = fft2(filt_tf_down_tuned);
         fbank_sr_domain(i,num_rate_ctrs+j,:,:) = filt_sr_down;
                         
    end
end

end









