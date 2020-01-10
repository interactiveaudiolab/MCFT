%% Example(2):
 % In this example, the Multi-resolution Common Fate Transform (MCFT)
 % of an audio signal is computed, and then the 4d-representation is
 % converted back to time domain. 
 
addpath([cd(cd('..')),'/mcft']); 
 
 %% Load the audio signal (trombone sample from the Philharmonia Orchestra
 % dataset):
 
 [aud_sig,fs] = audioread('trombone_D4_tremolo.wav');
 sig_dur = length(aud_sig)/fs;
 
 %% Compute and plot the CQT of the audio signal:
 % (see the CQT toolbox at http://www.cs.tut.fi/sgn/arg/CQT)
 
fmin = 27.5*2^(16/12); %C2
fmax = 27.5*2^(76/12); %C7
fres = 24; % bins per octave
gamma = 0;

sig_cq_struct = cqt(aud_sig, fres, fs, fmin, fmax,'gamma',gamma,'rasterize','full');
sig_cqt = sig_cq_struct.c;
freq_vec = sig_cq_struct.fbas;
time_frame_vec = linspace(0,sig_dur,size(sig_cqt,2));

figure
imagesc(time_frame_vec,freq_vec,20*log10(abs(sig_cqt)))
set(gca,'ydir','normal','fontsize',14,'fontweight','bold');
ylabel('\bf frequency (Hz)')
xlabel('\bf time (sec)')
colormap(parula)
box off
title('Constant-Q Treansform of the audio signal')

%% Compute the MCFT of the audio signal

cqt_params_in = struct('fs',fs,'fmin',fmin,'fmax',fmax,'fres',fres,'gamma',gamma);

[mcft_out,cqt_params_out,fbank_sr_domain] = mcft(aud_sig,cqt_params_in);

%% Reconstruct the time-domain signal

est_aud_sig = inv_mcft(mcft_out,cqt_params_out,fbank_sr_domain);

%% Compare the reconstructed signal to the original signal

% Error to signal ratio
reconst_err = 20*log10(norm(est_aud_sig - aud_sig)/norm(est_aud_sig));
disp(['Error to signal ratio = ' num2str(reconst_err) ' dB']);

% plot the original and reconstructed signals
time_vec=(0:length(aud_sig)-1)/fs;

figure
subplot(211)
plot(time_vec,aud_sig)
axis tight
ylabel('\bf x(t)')
xlabel('\bf time (sec)')
title('\bf original signal')

subplot(212)
plot(time_vec,est_aud_sig)
axis tight
ylabel('\bf xrec(t)')
xlabel('\bf time (sec)')
title('\bf reconstructed signal')










