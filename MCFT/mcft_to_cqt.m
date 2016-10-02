function X_hat=mcft_to_cqt(Z,H)

% This function reconstructs the time-frequency representation (CQT) of 
% an audio signal through inverse filtering given the 4-dimensioanl
% MCFT representation and the scale-rate domain filterbank.
% 
% Inputs: 
% Z: 4d matrix containing the MCFT coefficients
% H: 4d matrix containing a bank of scale-rate filters
%
% Output:
% X_hat: 2d matrix containing the reconstructed time-frequency
%        representation 
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)


%% Dimensions

[Ns,Nr,nfft_s,nfft_r]=size(H);

%% MCFT to TF representation

Hsr_sum=0;
Xsr_sum=0;
for i=1:Ns
    for j=1:Nr
         
         Z_temp=squeeze(Z(i,j,:,:));
         Xsr_temp=fft2(Z_temp,nfft_s,nfft_r);     
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
