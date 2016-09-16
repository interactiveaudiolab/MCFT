function  mcft_out=cqt_to_mcft(X,H) 

% This function receives the frequency-time representation (CQT) of 
% an audio signal (complex in general) and generates the 4-dimensional 
% representation (scale,rate, frequency, time) by 2d filtering based 
% on the cortical part of Chi's auditory model. 
%
% Inputs:
% X: 2d matrix contatining the (complex) t-f representation of an 
%    audio signal (log scale frequency, e.g. CQT)
% H: 4d matrix containing a bank of scale-rate filters
%
% Ouput: 
% Z: 4d matrix containing MCFT coefficients
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Dimensions

[Nf,Nt]=size(X);
Nf=Nf+mod(Nf,2); % next even number
Nt=Nt+mod(Nt,2);

[Ns,Nr,nfft_s,nfft_r]=size(H);

%% FT representation to MCFT:

% 2D-Fourier transform of the TF representation
Xsr=fft2(X,nfft_s,nfft_r);

mcft_out = zeros(Ns, Nr, Nf, Nt);
for i = 1:Ns
   for j=1:Nr
      Hsr_temp=squeeze(H(i,j,:,:));
      XHsr=Xsr.*Hsr_temp;   % filter the signal in scale-rate domain
      XHft=ifft2(XHsr); % convert back to the frequency-time domain
      mcft_out(i,j,:,:)=XHft(1:Nf,1:Nt);
           
   end
end

end


