%% Example(1):
 % In this example, we generate and plot upward and downward filter
 % impulse responses (STRFs) with a scale of 1 cycle per octave and
 % rate of 4 Hz.

addpath([cd(cd('..')),'/MCFT']); 
   
S = 1; R = 4;

% parameters
SRF=24; % spectral ripple frequency, or, # of bins per octave
Ls=3*SRF; % # of frequency channels of the spectrogram (3 octaves)
FPS= 150; % # of frames per second
Lr=FPS; % # of time frames of the spectrogram (1 second)
beta=3.5; % time constant
 
S_params=struct('hslen',Ls,'samprate_spec',SRF,'type','bandpass');
R_params=struct('time_const',beta,'hrlen',Lr,'samprate_temp',FPS,'type','bandpass');

start = tic;
[h_up,~]=gen_hsr(S,R,S_params,R_params,'up');
[h_down,~]=gen_hsr(S,R,S_params,R_params,'down');
stop = toc(start);
disp(['computation time: ',num2str(stop)]);

h_up=fftshift(h_up,1);
h_down=fftshift(h_down,1);

% frequency and time vectors for plotting
w = (0:Ls-1)'/SRF; % in number of octaves
t = (0:Lr-1)'/FPS; % in seconds

% plots
figure;
subplot(211)
imagesc(t,w,real(h_up))
set(gca,'ydir','normal','fontsize',14,'fontweight','bold');
ylabel('\bf octave number')
xlabel('\bf time (sec)')
colormap(parula)
box off
title('upward ripple')

subplot(212)
imagesc(t,w,real(h_down))
set(gca,'ydir','normal','fontsize',14,'fontweight','bold');
ylabel('\bf octave Number')
xlabel('\bf time (sec)')
colormap(parula)
box off
title('downward ripple')


%% test the filterbank

% generate a complex signal for phase modulation
fs = 8000;
t = (0:2*fs-1)/fs;

fmin = 27.5*2^(0/12);
fmax = 27.5*2^(87/12);
fres = 24; % bins per octave
gamma = 0;

signal = cos(2*pi*220*t)+cos(2*pi*440*t)+cos(2*pi*880*t);
Xcq = cqt(signal, fres, fs, fmin, fmax, 'gamma',gamma);
comp_specgram = Xcq.c;

% filterbank parameters
scale_ctrs = 2.^(-2:3);
rate_ctrs = 2.^(-2:0.125:6);
[nfft_s,nfft_r] = size(comp_specgram);

params = struct('samprate_spec',fres,'samprate_temp',nfft_r,'time_const',2);

start = tic;
[h_out,H_out] = gen_fbank_hsr(scale_ctrs,rate_ctrs,nfft_s,nfft_r,params); %,comp_specgram);
stop = toc(start);
disp(['computation time:',num2str(stop)]);


%% test cqt_to_mcft

start = tic;
mcft_out = cqt_to_mcft(comp_specgram,H_out);
stop = toc(start);
disp(['computation time:',num2str(stop)]);


%% test mcft

del_cqt_phase = 0;
cqt_params_in = struct('fs',fs,'fmin',fmin,'fmax',fmax,...
    'fres',fres,'gamma',gamma);
start = tic;
[mcft_out,cqt_params_out,H]=mcft(signal,cqt_params_in,del_cqt_phase);
stop = toc(start);
disp(['computation time:',num2str(stop)]);

%% test mcft to cqt

start = tic;
X_hat=mcft_to_cqt(mcft_out,H);
stop = toc(start);
disp(['computation time:',num2str(stop)]);

%% test inv_imcft

start = tic;
x_hat=inv_mcft(mcft_out,cqt_params_out,H);
stop = toc(start);
disp(['computation time:',num2str(stop)]);

rec_err = 20*log10(norm(x_hat.'-signal)/norm(signal));





