%% Example(1):
 % In this example, upward and downward filter impulse responses (STRFs)
 % with a scale of 1 cycle per octave and a rate of 4 cycles per second
 % are generated and plotted. 

addpath([cd(cd('..')),'/mcft']); 
   
scale_ctr = 1; 
rate_ctr = 4;

% parameters
spec_samprate = 24; % spectral ripple frequency, or, # of bins per octave
spec_nfft = 3 * spec_samprate; % # of frequency channels of the spectrogram (3 octaves)
temp_samprate = 150; % # of frames per second
temp_nfft = temp_samprate; % # of time frames of the spectrogram (1 second)
beta = 3.5; % time constant
 
scale_params = struct('spec_filt_len',spec_nfft,'spec_samprate',spec_samprate,...
    'spec_filt_type','bandpass');
rate_params = struct('temp_filt_len',temp_nfft,'temp_samprate',temp_samprate,...
    'temp_filt_type','bandpass','time_const',beta);

start = tic;
[filt_tf_up,~] = gen_filt_scale_rate(scale_ctr,rate_ctr,scale_params,...
    rate_params,'up');
[filt_tf_down,~] = gen_filt_scale_rate(scale_ctr,rate_ctr,scale_params,...
    rate_params,'down');
stop = toc(start);
disp(['computation time: ',num2str(stop)]);

filt_tf_up = fftshift(filt_tf_up,1);
filt_tf_down = fftshift(filt_tf_down,1);

% frequency and time vectors for plotting
freq_vec = (0:spec_nfft - 1)'/spec_samprate; % in number of octaves
time_frame_vec = (0:temp_nfft - 1)'/temp_samprate; % in seconds

% plots
figure
subplot(211)
imagesc(time_frame_vec,freq_vec,real(filt_tf_up))
set(gca,'ydir','normal','fontsize',14,'fontweight','bold');
ylabel('\bf octave number')
xlabel('\bf time (sec)')
colormap(parula)
box off
title('upward ripple')

subplot(212)
imagesc(time_frame_vec,freq_vec,real(filt_tf_down))
set(gca,'ydir','normal','fontsize',14,'fontweight','bold');
ylabel('\bf octave number')
xlabel('\bf time (sec)')
colormap(parula)
box off
title('downward ripple')


%% test the filterbank

% generate a complex signal for phase modulation
fs = 8000;
time_vec = (0:2*fs-1)/fs;

fmin = 27.5*2^(0/12);
fmax = 27.5*2^(87/12);
fres = 24; % bins per octave
gamma = 0;

signal = cos(2 * pi * 220 * time_vec) + cos(2 * pi * 440 * time_vec) + ...
    cos(2 * pi * 880 * time_vec);
sig_cq_struct = cqt(signal, fres, fs, fmin, fmax, 'gamma',gamma);
complex_specgram = sig_cq_struct.c;

% filterbank parameters
scale_ctrs = 2.^(1:5);
rate_ctrs = 2.^(-2:0.5:3);
[nfft_scale,nfft_rate] = size(complex_specgram);

scale_filt_params = struct('scale_ctrs',scale_ctrs,'nfft_scale',nfft_scale,...
    'spec_samprate',spec_samprate);
rate_filt_params = struct('rate_ctrs',rate_ctrs,'nfft_rate',nfft_rate,...
    'temp_samprate',temp_samprate,'time_const',beta);

start = tic;
[fbank_tf_domain,fbank_sr_domain] = gen_fbank_scale_rate(scale_filt_params,...
    rate_filt_params,sig_cqt);
stop = toc(start);
disp(['computation time:',num2str(stop)]);


%% test cqt_to_mcft

start = tic;
mcft_out = cqt_to_mcft(complex_specgram,fbank_sr_domain);
stop = toc(start);
disp(['computation time:',num2str(stop)]);


%% test mcft

del_cqt_phase = 0;
cqt_params_in = struct('fs',fs,'fmin',fmin,'fmax',fmax,...
    'fres',fres,'gamma',gamma);
start = tic;
[mcft_out,cqt_params_out,fbank_sr_domain] = mcft(signal,cqt_params_in,...
    'del_cqt_phase',del_cqt_phase,'predef_scale_ctrs',scale_ctrs,...
    'predef_rate_ctrs',rate_ctrs,'time_const',beta);
stop = toc(start);
disp(['computation time:',num2str(stop)]);


%% test mcft to cqt

start = tic;
est_cqt = mcft_to_cqt(mcft_out,fbank_sr_domain);
stop = toc(start);
disp(['computation time:',num2str(stop)]);



%% test inv_imcft

start = tic;
est_sig = inv_mcft(mcft_out,cqt_params_out,fbank_sr_domain);
stop = toc(start);
disp(['computation time:',num2str(stop)]);

rec_err = 20*log10(norm(est_sig.' - signal)/norm(signal));
disp(['reconstruction error:',num2str(rec_err)]);




