function H_out=gen_fbank_hsr_old(sv,rv,nfft_s,nfft_r,params,X)

% This function generates a scale-rate domain (analytic) filterbank. 
% The filterbank will be tuned to the passband of a target signal if
% specified.
% 
% Inputs: 
% sv: vector containing a range of scale values
% rv: vector containing a range of rate values (only +)
% nfft_s: number of scale fft points 
% nfft_r: number of rate fft points
% params: structure array containing filter parameters including:
%        1.ripple_freq 2.time_const 3.frame_per_sec 
% X: (optional) Nf by Nt matrix containing a 2d complex signal. 
%     If X is not provided the function will return the original
%     set of filters.   
% 
% Output: 
% H_out: Ns*Nr*(2Nf)*(2Nt) matrix containing the filterbank
% 
% Note: the first and last filters in s and r ranges are assumed
%       to be lowpass and highpass respectively

%% Input check
tune_filter=0;
if nargin==6
    tune_filter=1;
end

%% Parameters and dimensions

rv_full=[-rv(end:-1:1),rv]; % full rate vector (-/+)

% Nf=2^nextpow2(Nf); 
% Nt=2^nextpow2(Nt);
nfft_s=nfft_s+mod(nfft_s,2);
nfft_r=nfft_r+mod(nfft_r,2);

Ns=length(sv);
Nr=length(rv_full);

beta=params.time_const;
SRF=params.ripple_freq;
FPS=params.frame_per_sec;

s_params=struct('hslen',nfft_s,'ripple_freq',SRF);
r_params=struct('time_const',beta,'hrlen',nfft_r,'frame_per_sec',FPS);


%% Generate the filterbank

% Filter tuning factor (can be thought as a pre-filtering stage)
if tune_filter
   Xft=ifft2(fft2(X,nfft_s,nfft_r)); % for dimension adjustment
   %Xft_mask=abs(Xft)/max(abs(Xft(:)));
   %mod_ang=angle(Xft).*Xft_mask;
   %mod_ang=mod_ang*pi/max(abs(mod_ang(:)));
   h_factor=exp(1j.*angle(Xft));
    
%   H_factor=fft2(abs(X),nfft_s,nfft_r)./(fft2(X,nfft_s,nfft_r)+eps);
%   H_factor=(abs(H_factor)/max(H_factor(:))).*exp(1j*angle(H_factor));
else
   h_factor=1; 
    
%   H_factor=1;
end

H_out = zeros(Ns, Nr, nfft_s, nfft_r);
for i = 1:Ns
    %i
    
    if sv(i)==sv(1) 
        s_params.type='lowpass'; 
        %s_params.type='bandpass'; 

    elseif sv(i)==sv(end)
        s_params.type='highpass'; 
        %s_params.type='bandpass'; 

    else
        s_params.type='bandpass';
    end
    
    for j=1:Nr
    
         if abs(rv_full(j))==rv(1)
           r_params.type='lowpass';
           %r_params.type='bandpass';
         elseif abs(rv_full(j))==rv(end)
           r_params.type='highpass';
           %r_params.type='bandpass';

         else
           r_params.type='bandpass';
         end
         
         % consider the sign of rate values
         if sign(rv_full(j))<0
             r_params.r_support='neg';
         elseif sign(rv_full(j))>0
             r_params.r_support='pos';
         end
         
         %[~,Hsr]=gen_hsr(sv(i),abs(rv_full(j)),s_params,r_params);
         [hsr,~]=gen_hsr_old(sv(i),abs(rv_full(j)),s_params,r_params);
         hsr_tuned=hsr.*h_factor;
         Hsr_temp=fft2(hsr_tuned);

         %Hsr_temp=Hsr.*H_factor;   % tune the filter 
         H_out(i,j,:,:)=Hsr_temp;
        
    end
end










