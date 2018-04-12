%% Example(2):
 % In this example, we compute the Multi-resolution Common Fate Transform
 % (MCFT) of an audio signal, and then convert the 4d-representation back
 % to time domain. 
 
addpath([cd(cd('..')),'/MCFT']); 
 
 %% Load the audio signal (trombone sample from the Philharmonia Orchestra
 % dataset):
 
 [x,fs]=audioread('trombone_D4_tremolo.wav');
 dur_x=length(x)/fs;
 
 %% Compute and plot the CQT of the audio signal:
 % (see the CQT toolbox at http://www.cs.tut.fi/sgn/arg/CQT)
 
fmin = 27.5*2^(16/12); %C2
fmax = 27.5*2^(76/12); %C7
fres = 24; % bins per octave

Xcq = cqt(x, fres, fs, fmin, fmax,'gamma',0,'rasterize','full');
X=Xcq.c;
fvec=Xcq.fbas;
tvec=linspace(0,dur_x,size(X,2));

figure
imagesc(tvec,fvec,20*log10(abs(X)))
set(gca,'ydir','normal','fontsize',14,'fontweight','bold');
ylabel('\bf frequency (Hz)')
xlabel('\bf time (sec)')
colormap(parula)
box off
title('Constant-Q Treansform of the audio signal')

%% Compute the MCFT of the audio signal

cqt_params_in=struct('fs',fs,'fmin',fmin,'fmax',fmax,'fres',fres);

[Z,cqt_params_out,H]=mcft(x,cqt_params_in);

%% Reconstruct the time-domain signal

x_hat=inv_mcft(Z,cqt_params_out,H);

%% Compare the reconstructed signal to the original signal

% Error to signal ratio
Err=20*log10(norm(x_hat-x)/norm(x_hat));
disp(['Error to signal ratio = ' num2str(Err) ' dB']);

% plot the original and reconstructed signals
tvec=(0:length(x)-1)/fs;
figure
subplot(211)
plot(tvec,x)
axis tight
ylabel('\bf x(t)')
xlabel('\bf time (sec)')
title('\bf original signal')

subplot(212)
plot(tvec,x_hat)
axis tight
ylabel('\bf xrec(t)')
xlabel('\bf time (sec)')
title('\bf reconstructed signal')










