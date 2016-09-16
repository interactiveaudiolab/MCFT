%% Example(1):
 % In this example, we generate and plot upward and downward filter
 % impulse responses (STRFs) with a scale of 1 cycle per octave and
 % rate of 4 Hz.

addpath([cd(cd('..')),'/mcft']); 
   
S=1; R=4;

% parameters
SRF=24; % spectral ripple frequency, or, # of bins per octave
Ls=3*SRF; % # of frequency channels of the spectrogram (3 octaves)
FPS=150; % # of frames per second
Lr=FPS; % # of time frames of the spectrogram (1 second)
beta=3.5; % time constant
 
S_params=struct('hslen',Ls,'ripple_freq',SRF,'type','bandpass');
R_params=struct('time_const',beta,'hrlen',Lr,'frame_per_sec',FPS,'type','bandpass');

[h_up,~]=gen_hsr(S,R,S_params,R_params,'up');
[h_down,~]=gen_hsr(S,R,S_params,R_params,'down');

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

