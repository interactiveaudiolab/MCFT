function [h_out,H_out] = gen_fbank_hsr(scale_ctrs,rate_ctrs,nfft_s,nfft_r,params,comp_specgram)

% This function generates a scale-rate domain bank of up-/down-ward,
% filters. The filterbank will be tuned to the passband of a target 
% signal if specified.
% 
% Inputs: 
% scale_ctrs: vector containing filter centers along the scale axis
% rate_ctrs: vector containing filter centers along the rate axis
% nfft_s: number of scale fft points 
% nfft_r: number of rate fft points
% params: structure array containing filter parameters including:
%        1.samprate_spec 2.samprate_temp 3.time_const 
% comp_specgram: (optional) Nf by Nt matrix containing a 2d complex spectrogram. 
%     If provided, the function will return a filter bank that is modulated with 
%     the phase of the spectrogram, otherwise, the function will return the 
%     original set of filters.   
% 
% Output: 
% h_out: Ns*(2*Nr)*nfft_s*nfft_r matix containing the time-freq-domain filter bank
% H_out: Ns*(2*Nr)*nfft_s*nfft_r matrix containing the scale-rate-domain filter bank
% note: nfft_s>=Nf and nfft_r>=Nt where [Nf,Nt]=size(X)
%
% Note: the first and last filters in scale and rate ranges are assumed
%       to be lowpass and highpass respectively
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)


%% Input check
tune_filter = 0;
if nargin == 6
    tune_filter = 1;
end

%% Parameters and dimensions

nfft_s = nfft_s+mod(nfft_s,2); % set nfft to the next even number
nfft_r = nfft_r+mod(nfft_r,2);

num_scale_ctrs = length(scale_ctrs);
num_rate_ctrs = length(rate_ctrs);

beta = params.time_const;
samprate_spec = params.samprate_spec;
samprate_temp = params.samprate_temp;

scale_params=struct('hslen',nfft_s,'samprate_spec',samprate_spec);
rate_params=struct('time_const',beta,'hrlen',nfft_r,'samprate_temp',samprate_temp);


%% Generate the filterbank

% Filter tuning factor (can be thought as a pre-filtering stage)
if tune_filter
   Xft=ifft2(fft2(comp_specgram,nfft_s,nfft_r)); % for dimension adjustment
%    [Nf,Nt]=size(comp_specgram);
%    Xft=zeros(nfft_s,nfft_r);
%    Xft(1:Nf,1:Nt)=comp_specgram;
   h_factor=exp(1j.*angle(Xft));
else
   h_factor=1; 
end

h_out = zeros(num_scale_ctrs, 2*num_rate_ctrs, nfft_s, nfft_r);
H_out = zeros(num_scale_ctrs, 2*num_rate_ctrs, nfft_s, nfft_r);
for i = 1:num_scale_ctrs
    
    if scale_ctrs(i) == scale_ctrs(1) 
        scale_params.type = 'lowpass'; 
        
    elseif scale_ctrs(i) == scale_ctrs(end)
        scale_params.type = 'highpass'; 
    else
        scale_params.type = 'bandpass';
    end
    
    for j = 1:num_rate_ctrs
    
         if rate_ctrs(j) == rate_ctrs(1)
           rate_params.type = 'lowpass';
         elseif rate_ctrs(j) == rate_ctrs(end)
           rate_params.type = 'highpass';
         else
           rate_params.type = 'bandpass';
         end
         
         % generate two analytic filters (one upward and one downward)
         % for the current (S,R) values
         
         % upward
         [hsr_up,~] = gen_hsr(scale_ctrs(i),rate_ctrs(j),scale_params,rate_params,'up');
         hsr_up_tuned = hsr_up.*h_factor;
         h_out(i,num_rate_ctrs-j+1,:,:) = hsr_up_tuned;
         Hsr_up = fft2(hsr_up_tuned);
         H_out(i,num_rate_ctrs-j+1,:,:) = Hsr_up;
         
         % downward
         [hsr_down,~] = gen_hsr(scale_ctrs(i),rate_ctrs(j),scale_params,rate_params,'down');
         hsr_down_tuned = hsr_down.*h_factor;
         h_out(i,num_rate_ctrs+j,:,:) = hsr_down_tuned;
         Hsr_down = fft2(hsr_down_tuned);
         H_out(i,num_rate_ctrs+j,:,:) = Hsr_down;
                         
    end
end

end









