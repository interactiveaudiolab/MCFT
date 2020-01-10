%% Test the gen_fbank_hsr function

clc; clear; close all

%% signal

load('unison_mix_D4.mat','all_sources','all_mixtures','all_names');

num_samp=size(all_sources,1);

example_num = 1;

s1 = all_sources{example_num}(1,:);
s2 = all_sources{example_num}(2,:);
x = all_mixtures{example_num};

%% CQT

% cqt parameters
fs=22050;
fmin = 27.5*2^((15)/12); %C2
fres = 24;
gamma = 0; 
fmax = 27.5*2^((75)/12); %C7

% sources 
Scq1 = cqt(s1, fres, fs, fmin, fmax,'gamma',gamma,'rasterize','full');
S1_cqt=Scq1.c;

Scq2 = cqt(s2, fres, fs, fmin, fmax,'gamma',gamma,'rasterize','full');
S2_cqt=Scq2.c;

% mixture
Xcq = cqt(x, fres, fs, fmin, fmax,'gamma',gamma,'rasterize','full');
X_cqt=Xcq.c;

[Nf,Nt]=size(X_cqt);
freq_vec=Xcq.fbas;
dur_x = length(x)/fs;
time_vec=linspace(0,dur_x,Nt);

fmat=repmat(freq_vec,1,Nt);
Xph=unwrap(angle(X_cqt),[],2)./fmat;
Xcqt_fnorm=abs(X_cqt).*exp(1j*Xph);

%% common parameters

SRF=fres; 
FPS=floor(Nt/dur_x);

nfft_s = Nf;
nfft_r = Nt;

smax=nextpow2(SRF/2)-1;
rmax=nextpow2(FPS/2)-1;

sv = 2.^[(-3:1:smax),3.3]; 
rv = 2.^[(-3:0.5:rmax),5.75]; 

beta = 1;

%% parameters of gen_fbank_hsr_old

H_params_old=struct('time_const',beta,'ripple_freq',SRF,'frame_per_sec',FPS);

%% parameters of gen_fbank_hsr

H_params=struct('samprate_spec',SRF,'samprate_temp',FPS,'time_const',beta);

%% compare filter banks

H_old=gen_fbank_hsr_old(sv,rv,nfft_s,nfft_r,H_params_old,X_cqt); 

[~,H]=gen_fbank_hsr(sv,rv,nfft_s,nfft_r,H_params,X_cqt); 

disp(['total error = ', num2str(sum(abs(H_old(:)-H(:))))])

%% compute the filter bank through mcft function

cqt_params_in=struct('fs',fs,'fmin',fmin,'fmax',fmax,'fres',fres);
filt2d_params = struct('scale_ctrs',sv,'rate_ctrs',rv,'time_const',beta);

[X_mcft,~,H_mcft]=mcft(x,cqt_params_in,filt2d_params); 

disp(['total error = ', num2str(sum(abs(H_old(:)-H_mcft(:))))])






