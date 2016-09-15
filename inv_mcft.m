function x_hat=inv_mcft(Z,cqt_params,H)

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
% Z: 4d matrix containing MCFT coefficients
% H: 4d matrix containing the scale-rate filter bank
%
% Ouput:
% x_hat: vector containing the reconstructed time-domain signal
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwesterin.edu)

%% MCFT to CQT

X_hat=mcft_to_cqt(Z,H); % reconstructed CQT of the signal

%% CQT to time-domain signal

Xcq=cqt_params;
Xcq.c=X_hat;
x_hat=icqt(Xcq);


end