function  mcft_out=cqt_to_mcft(Xcqt,H) 

% This function receives the frequency-time representation (CQT) of 
% an audio signal (complex in general) and generates the 4-dimensional 
% representation (scale,rate, frequency, time) by 2d filtering based 
% on the cortical part of Chi's auditory model. 
%
% Inputs:
% Xcqt: 2d matrix contatining the (complex) t-f representation of an 
%    audio signal (log scale frequency, e.g. CQT)
% H: 4d matrix containing a bank of scale-rate filters
%
% Ouput: 
% mcft_out: 4d matrix containing MCFT coefficients
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Dimensions

[Nf,Nt] = size(Xcqt);
Nf = Nf+mod(Nf,2); % next even number
Nt = Nt+mod(Nt,2);

[num_scale_ctrs,num_rate_ctrs,nfft_s,nfft_r] = size(H);

%% FT representation to MCFT:

% 2D-Fourier transform of the TF representation
Xsr = fft2(Xcqt,nfft_s,nfft_r);

mcft_out = zeros(num_scale_ctrs, num_rate_ctrs, nfft_s, nfft_r);
for i = 1:num_scale_ctrs
   for j = 1:num_rate_ctrs
      Hsr_temp = squeeze(H(i,j,:,:));
      XHsr = Xsr.*Hsr_temp;   % filter the signal in scale-rate domain
      XHft = ifft2(XHsr); % convert back to the frequency-time domain
      mcft_out(i,j,:,:)=XHft;%(1:Nf,1:Nt); % remove the zero padding
           
   end
end

end


