function H_out=gen_fbank_hsr(sv,rv,nfft_s,nfft_r,params,X)

% This function generates a scale-rate domain bank of up-/down-ward,
% filters. The filterbank will be tuned to the passband of a target 
% signal if specified.
% 
% Inputs: 
% sv: vector containing a range of scale values
% rv: vector containing a range of rate values 
% nfft_s: number of scale fft points 
% nfft_r: number of rate fft points
% params: structure array containing filter parameters including:
%        1.ripple_freq 2.frame_per_sec 3.time_const 
% X: (optional) Nf by Nt matrix containing a 2d complex signal. 
%     If X is not provided, the function will return the original
%     set of filters.   
% 
% Output: 
% H_out: Ns*(2*Nr)*Nf*Nt matrix containing the filterbank
% 
% Note: the first and last filters in s and r ranges are assumed
%       to be lowpass and highpass respectively
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)


%% Input check
tune_filter=0;
if nargin==6
    tune_filter=1;
end

%% Parameters and dimensions

nfft_s=nfft_s+mod(nfft_s,2); % set nfft to the next even number
nfft_r=nfft_r+mod(nfft_r,2);

Ns=length(sv);
Nr=length(rv);

beta=params.time_const;
SRF=params.ripple_freq;
FPS=params.frame_per_sec;

S_params=struct('hslen',nfft_s,'ripple_freq',SRF);
R_params=struct('time_const',beta,'hrlen',nfft_r,'frame_per_sec',FPS);


%% Generate the filterbank

% Filter tuning factor (can be thought as a pre-filtering stage)
if tune_filter
   Xft=ifft2(fft2(X,nfft_s,nfft_r)); % for dimension adjustment
   h_factor=exp(1j.*angle(Xft));
else
   h_factor=1; 
end

H_out = zeros(Ns, 2*Nr, nfft_s, nfft_r);
for i = 1:Ns
    
    if sv(i)==sv(1) 
        S_params.type='lowpass'; 
        
    elseif sv(i)==sv(end)
        S_params.type='highpass'; 
    else
        S_params.type='bandpass';
    end
    
    for j=1:Nr
    
         if rv(j)==rv(1)
           R_params.type='lowpass';
         elseif rv(j)==rv(end)
           R_params.type='highpass';
         else
           R_params.type='bandpass';
         end
         
         % generate two analytic filters (one upward and one downward)
         % for the current (S,R) values
         
         % upward
         [hsr_up,~]=gen_hsr(sv(i),rv(j),S_params,R_params,'up');
         hsr_up_tuned=hsr_up.*h_factor;
         Hsr_up=fft2(hsr_up_tuned);
         H_out(i,Nr-j+1,:,:)=Hsr_up;
         
         % downward
         [hsr_down,~]=gen_hsr(sv(i),rv(j),S_params,R_params,'down');
         hsr_down_tuned=hsr_down.*h_factor;
         Hsr_down=fft2(hsr_down_tuned);
         H_out(i,Nr+j,:,:)=Hsr_down;
         
         %rv(j)
         %Nr-j+1
         %Nr+j
                 
    end
end

end









