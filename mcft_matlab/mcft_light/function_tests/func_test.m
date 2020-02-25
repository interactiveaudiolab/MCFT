%% In this script toolbox funtions of the mcft_light toolbox are tested
clc; clear; close all

%% generate a test signal

% generate the signal   
fs = 2*4000;
fund_freq = 27.5*2^((40-1)/12); 

time_vec = 0:1/fs:1-1/fs; 
num_harmonics = 5;   
harmonic_freq = fund_freq * (1:num_harmonics)';
modulation_amp = 20*(1:num_harmonics)'; 
modulation_freq = 3;
sig_1 = cos(2 * pi * fund_freq * time_vec);
sig_2 = cos(2 * pi * fund_freq * time_vec + ...
    modulation_amp(1) * cos(2 * pi * modulation_freq * time_vec));
sig_harmonic_1 = sum(cos(2 * pi * harmonic_freq * time_vec)); 
sig_harmonic_2 = sum(cos(2 * pi * harmonic_freq * time_vec + ...
    modulation_amp * cos(2 * pi * modulation_freq * time_vec)));
sig = sig_harmonic_1 + sig_harmonic_2;

sig_dur = length(sig)/fs;

% cqt of the signal 
fmin = 27.5*2^(0/12); %C2
fmax = 27.5*2^(87/12); %C8
fres = 48; % bins per octave

sig_cq_struct = cqt(sig, fres, fs, fmin, fmax,'gamma',50,'rasterize','full');
sig_cqt = sig_cq_struct.c;

[n_freq,n_time_frame] = size(sig_cqt);
freq_vec = sig_cq_struct.fbas;
time_frame_vec = linspace(0,sig_dur,n_time_frame);

fvec_khz = freq_vec/1000;
win2d = 1; %blackman(n_freq) * blackman(n_time_frame).';

scale_nfft = n_freq;
rate_nfft = n_time_frame;
sig_sr_domain = fftshift(fft2(abs(win2d .* sig_cqt),scale_nfft,rate_nfft));

% colormap
ngrays=400;
lgrays=zeros(ngrays,3);
for i=1:ngrays
    lgrays(i,:) = 1-(i/ngrays)^0.5;
end 

figure(1)
subplot(211)
surf(time_frame_vec,freq_vec,abs(sig_cqt),'Linestyle','none')
view(0,90)
set(gca,'yscale','log');
axis tight
yticks([250,500,1000,2000])
yticklabels({'250','500','1000','2000'})
ylab=ylabel('Frequency (Hz)');
set(ylab, 'Units', 'Normalized', 'position',[-0.1, 0.5, 1]);
xlabel('Time (sec)')
title('(a)')
colormap(lgrays)
colorbar off

subplot(212)
mesh(abs(sig_sr_domain))
set(gca,'ydir','normal')
axis tight

%% function gen_1d_fbank_light

% scale filterbank
spec_samprate = fres;
scale_ctr_min = 1;
scale_ctr_max = 2^(nextpow2(spec_samprate/2)-1);
scale_filt_res = 1;

[scale_fbank,scale_filt_ctrs,scale_ctr_posit,scale_filt_len] = ...
    gen_1d_fbank_light('gabor_fourier',scale_ctr_min,scale_ctr_max,...
    scale_filt_res,spec_samprate,scale_nfft);

n_scale_ctr = length(scale_fbank);
scale_hpass_idx = floor(n_scale_ctr/2) + 1;

figure(2)
subplot(scale_hpass_idx,1,1)
stem(scale_fbank{1})
grid on; axis tight
title('lowpass scale filter')

subplot(scale_hpass_idx,1,scale_hpass_idx)
stem(scale_fbank{scale_hpass_idx})
grid on; axis tight
title('highpass scale filter')

for i = 2:scale_hpass_idx-1
    subplot(scale_hpass_idx,1,i)
    stem(scale_fbank{i})
    grid on; axis tight
    title('bandpass scale filter')
    
end

% rate filterbank
temp_samprate = n_time_frame/sig_dur;
rate_ctr_min = 1;
rate_ctr_max = 2^(nextpow2(temp_samprate/2)-1);
rate_filt_res = 1;

[rate_fbank,rate_filt_ctrs,rate_ctr_posit,rate_filt_len] = ...
    gen_1d_fbank_light('gammatone_fourier',rate_ctr_min,rate_ctr_max,...
    rate_filt_res,temp_samprate,rate_nfft);

n_rate_ctr = length(rate_fbank);
rate_hpass_idx = floor(n_rate_ctr/2) + 1;

figure(3)
subplot(rate_hpass_idx,1,1)
stem(abs(rate_fbank{1}))
grid on; axis tight
title('lowpass rate filter')

subplot(rate_hpass_idx,1,rate_hpass_idx)
stem(abs(rate_fbank{rate_hpass_idx}))
grid on; axis tight
title('highpass rate filter')

for i = 2:rate_hpass_idx-1
    subplot(rate_hpass_idx,1,i)
    stem(abs(rate_fbank{i}))
    grid on; axis tight
    title('bandpass scale filter')
    
end

%% function gen_2d_fbank_light

ctr_min = [scale_ctr_min, rate_ctr_min];
ctr_max = [scale_ctr_max, rate_ctr_max];
filt_res = [scale_filt_res, rate_filt_res];
nfft = [scale_nfft, rate_nfft];
samprate = [spec_samprate, temp_samprate];

fbank_params = struct('ctr_min',ctr_min,'ctr_max',ctr_max,'filt_res',filt_res,...
    'nfft',nfft,'samprate',samprate);

[fbank,filt_ctr_posit,filt_range] = gen_2d_fbank_light(fbank_params,'zpadding','off');

figure(4)
for i = 1:n_scale_ctr
    for j = 1:n_rate_ctr
        subplot(n_scale_ctr,n_rate_ctr,(i-1)*n_rate_ctr+j)
        imagesc(fftshift(abs(fbank{i,j})))
        set(gca,'ydir','normal')
        axis tight
        
        if i==1 && j==1
            title('R\_low')
            ylabel('S\_low')
        elseif i==1 && j~=rate_hpass_idx && j~=rate_hpass_idx+1
            title('R\_band')
        elseif i==1 && (j==rate_hpass_idx || j==rate_hpass_idx+1)
            title('R\_high')
        elseif j==1 && i~=scale_hpass_idx && i~=scale_hpass_idx+1
            ylabel('S\_band')
        elseif j==1 && (i==scale_hpass_idx || i==scale_hpass_idx+1)
            ylabel('S\_high')
            
        end        
    end
end

%% function cqt_to_mcft_light

% compute the mcft 
sig_cqt_in = sig_cqt;
mcft_out = cqt_to_mcft_light(sig_cqt_in,fbank,filt_range,filt_ctr_posit,'zpadding','on');

figure(5)
for i = 1:n_scale_ctr
    for j = 1:n_rate_ctr
        subplot(n_scale_ctr,n_rate_ctr,(i-1)*n_rate_ctr+j)
        imagesc(abs(mcft_out{i,j}))
        set(gca,'ydir','normal')
        axis tight
        
        if i==1 && j==1
            title('R\_low')
            ylabel('S\_low')
        elseif i==1 && j~=rate_hpass_idx && j~=rate_hpass_idx+1
            title('R\_band')
        elseif i==1 && (j==rate_hpass_idx || j==rate_hpass_idx+1)
            title('R\_high')
        elseif j==1 && i~=scale_hpass_idx && i~=scale_hpass_idx+1
            ylabel('S\_band')
        elseif j==1 && (i==scale_hpass_idx || i==scale_hpass_idx+1)
            ylabel('S\_high')
            
        end        
    end
end

%% function mcft_to_cqt_light

mcft_in = mcft_out;
sig_cqt_size = [n_freq,n_time_frame];

est_sig_cqt = mcft_to_cqt_light(mcft_in,fbank,filt_ctr_posit,filt_range,...
    sig_cqt_size);

cqt_reconst_err = 20 * log10(norm(sig_cqt_in(:) - est_sig_cqt(:))/norm(sig_cqt_in(:)));

disp(['CQT reconstruction error = ' num2str(cqt_reconst_err) ' dB']);


%% function mcft_light

cqt_params_in = struct('fs',fs,'fmin',fmin,'fmax',fmax,'fres',fres,'gamma',0);

% fbank parameters
fbank_params = struct('ctr_min',ctr_min,'ctr_max',ctr_max,'filt_res',filt_res);

start = tic;
[mcft_out,inv_bundle] = mcft_light(sig,cqt_params_in,fbank_params,...
    'zpadding','off','output_format','cell','del_cqt_phase',0);
stop = toc(start)


%% function inv_mcft_light

mcft_in = mcft_out;

est_sig = inv_mcft_light(mcft_in,inv_bundle);

reconst_err = 20*log10(norm(est_sig(:) - sig(:)) / norm(sig(:)));
disp(['Error to signal ratio = ' num2str(reconst_err) ' dB']);







