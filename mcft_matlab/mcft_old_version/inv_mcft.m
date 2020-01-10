function x_hat=inv_mcft(mcft_in,cqt_params,H)

% This function reconstructs a time-domain signal given its 
% Multi-resolution Common Fate Transform (MCFT). 
% The intermediary time-frequency domain representation is the 
% Constant-Q Transform (CQT) of the audio signal, which is computed 
% using the invertible and optimized CQT implementation proposed by 
% Schorkhuber et al.:
%
% Toolbox webpage: 
% http://www.cs.tut.fi/sgn/arg/CQT/
% Reference: 
% Schörkhuber et al. "A Matlab toolbox for efficient perfect 
% reconstruction time-frequency transforms with log-frequency resolution."
% Audio Engineering Society Conference: 53rd International Conference: 
% Semantic Audio. Audio Engineering Society, 2014. 
%
% Inputs: 
% mcft_in: 4d matrix containing MCFT coefficients
% cqt_params: structure array containing cqt parameters including:
%             Nf: number of frequency bins
%             Nt: number of time frames
% H: 4d matrix containing the scale-rate filter bank
%
% Ouput:
% x_hat: vector containing the reconstructed time-domain signal
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% MCFT to CQT

disp('Reconstructing the CQT...');
X_hat=mcft_to_cqt(mcft_in,H); % reconstructed CQT of the signal


%% CQT to time-domain signal

disp('Reconstructing the time-domain signal...');

Nf=cqt_params.Nf;
Nt=cqt_params.Nt;
cqt_params=rmfield(cqt_params,'Nf');
cqt_params=rmfield(cqt_params,'Nt');

Xcq=cqt_params;

% see observation note in the mcft function
%%%% phase normalization method %%%%
% % reconstruct the phase
% freq_vec=Xcq.fbas;
% fmat=repmat(freq_vec,1,Nt);
% X_hat = X_hat(1:Nf,1:Nt);
% Xph=unwrap(angle(X_hat),[],2).*fmat;
%%% the following line is interesting and worth exploring
% X_hat=abs(X_hat).*exp(1j*angle(X_hat)./(abs(X_hat)+eps)); %Xph);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xcq.c=X_hat(1:Nf,1:Nt);
x_hat=icqt(Xcq);


end