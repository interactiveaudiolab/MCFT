function X_hat=inv_mcft(mcft_in,H)

% This function reconstructs the time-frequency representation of 
% an audio signal given the 4-dimensioanl MCFT representation and the 
% scale-rate domain filterbank.
% 
% Inputs: 
% mcft_in: 4-dimensional matrix containing the MCFT coefficients
% H: 4d matrix containing a set of scale-rate filters
%
% Output:
% X_hat: 2-dimensioanl matrix containing the reconstructed time-freq.
%        representation (log scale frequency, e.g. CQT)

%% Dimensions

[Ns,Nr,nfft_s,nfft_r]=size(H);

%% MCFT to TF representation

Hsr_sum=0;
Xsr_sum=0;
for i=1:Ns
    for j=1:Nr
         
         mcft_temp=squeeze(mcft_in(i,j,:,:));
         Xsr_temp=fft2(mcft_temp,nfft_s,nfft_r);     
         Hsr_hat=squeeze(H(i,j,:,:));
         
         Xhsr=Xsr_temp.*conj(Hsr_hat);
         
         Hsr_sum=Hsr_sum+Hsr_hat.*conj(Hsr_hat);
         Xsr_sum=Xsr_sum+Xhsr;
         
    end
end
                  
% compute the ratio Xsr_sum/Hsr_sum
Xsr_ratio=Xsr_sum./(Hsr_sum+eps);

% compute X_hat
X_hat=ifft2(Xsr_ratio);

end
