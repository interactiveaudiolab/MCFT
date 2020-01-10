function est_signal = inv_mcft_refactored(mcft_in,cqt_params,fbank_sr_domain)

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
% fbank_sr_domain: 4d matrix containing the scale-rate filter bank
%
% Ouput:
% est_signal: vector containing the reconstructed time-domain signal
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% MCFT to CQT

disp('Reconstructing the CQT...');
est_sig_cqt = mcft_to_cqt_refactored(mcft_in,fbank_sr_domain); % reconstructed CQT of the signal


%% CQT to time-domain signal

disp('Reconstructing the time-domain signal...');

n_freq = cqt_params.n_freq;
n_time = cqt_params.n_time;
cqt_params = rmfield(cqt_params,'n_freq');
cqt_params = rmfield(cqt_params,'n_time');

sig_cq_struct = cqt_params;

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

sig_cq_struct.c = est_sig_cqt(1:n_freq,1:n_time);
est_signal = icqt(sig_cq_struct);


end