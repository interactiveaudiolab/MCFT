%% Compare refactored code to the original code
clc; clear; close all;

%% filt_default_ceters 

% original function
nfft_s = 122;
nfft_r = 102;
samprate_spec = 24;
samprate_temp = 204;

[scale_ctrs_org,rate_ctrs_org] = filt_default_centers(nfft_s,nfft_r,samprate_spec,samprate_temp);

% refactored function
scale_res = 1;
rate_res = 1;

scale_params = struct('scale_res',scale_res,'scale_nfft',nfft_s,'samprate_spec',samprate_spec);
rate_params = struct('rate_res',rate_res,'rate_nfft',nfft_r,'samprate_temp',samprate_temp);

[scale_ctrs_ref1,rate_ctrs_ref1] = filt_default_centers_refactored1(scale_params,rate_params);


scale_params = struct('filt_type','scale','filt_res',scale_res,'filt_nfft',nfft_s,'samprate',samprate_spec);
rate_params = struct('filt_type','rate','filt_res',rate_res,'filt_nfft',nfft_r,'samprate',samprate_temp);

scale_ctrs_ref = filt_default_centers_refactored(scale_params);
rate_ctrs_ref = filt_default_centers_refactored(rate_params);

fprintf('Original scale centers: \n')
disp(scale_ctrs_org)
fprintf('refactored scale centers: \n')
disp(scale_ctrs_ref1)
fprintf('refactored2 scale centers: \n')
disp(scale_ctrs_ref)
fprintf('Original rate centers: \n')
disp(rate_ctrs_org)
fprintf('refactored rate centers: \n')
disp(rate_ctrs_ref1)
fprintf('refactored2 rate centers: \n')
disp(rate_ctrs_ref)

%% gen_hsr

scale_ctr = 4;
rate_ctr = 8;
filt_dir = 'down';

% original function
scale_params_org = struct('hslen',nfft_s,'samprate_spec',samprate_spec,...
    'type','bandpass');
rate_params_org = struct('hrlen',nfft_r,'samprate_temp',samprate_temp,...
    'type','bandpass','time_const',1);
[h_org,H_org] = gen_hsr(scale_ctr,rate_ctr,scale_params_org,rate_params_org,filt_dir);

% refactored function
scale_params_ref = struct('spec_filt_len',nfft_s,'spec_samprate',samprate_spec,...
    'spec_filt_type','bandpass');
rate_params_ref = struct('temp_filt_len',nfft_r,'temp_samprate',samprate_temp,...
    'temp_filt_type','bandpass','time_const',1);
[h_ref,H_ref] = gen_filt_scale_rate(scale_ctr,rate_ctr,scale_params_ref,...
    rate_params_ref,filt_dir);

h_err = norm(h_org(:) - h_ref(:));
H_err = norm(H_org(:) - H_ref(:));

%% gen_fbank_hsr

scale_ctrs = [1,2,4];
rate_ctrs = [1,2,4,8];

% original function
params = struct('samprate_spec',samprate_spec,'samprate_temp',samprate_temp,...
    'time_const',1);
[h_org,H_org] = gen_fbank_hsr(scale_ctrs,rate_ctrs,nfft_s,nfft_r,params);

% refactored function
scale_filt_params = struct('scale_ctrs',scale_ctrs,'nfft_scale',nfft_s,...
    'spec_samprate',samprate_spec);
rate_filt_params = struct('rate_ctrs',rate_ctrs,'nfft_rate',nfft_r,...
    'temp_samprate',samprate_temp,'time_const',1);
[h_ref,H_ref] = gen_fbank_scale_rate(scale_filt_params,rate_filt_params);

h_err = norm(h_org(:) - h_ref(:));
H_err = norm(H_org(:) - H_ref(:));

 
%% mcft

% cqt params
% CQT parameters
fs = 8000;
fmin = 110;
fmax = 110*2^5;
fres = 24; 
gamma = 0; 
cqt_params_in = struct('fs',fs,'fmin',fmin,'fmax',fmax,'fres',fres,'gamma',gamma);

del_cqt_phase = 0;

tt = (0:0.5*fs-1)/fs;
x = cos(2*pi*440*tt).';

% original function
[mcft_out_org,cqt_params_out_org,H_org] = mcft(x,cqt_params_in,del_cqt_phase);


% refactored function
[mcft_out_ref,cqt_params_out_ref,H_ref] = mcft_refactored(x,cqt_params_in,'del_cqt_phase',...
    del_cqt_phase);

% errors
mcft_err = norm(mcft_out_org(:) - mcft_out_ref(:))/norm(mcft_out_org(:))

figure(1)
subplot(221)
imagesc(squeeze(abs(mcft_out_org(end,8,:,:))))
set(gca,'ydir','normal')
subplot(222)
imagesc(squeeze(abs(fftshift(H_org(end,8,:,:)))))
set(gca,'ydir','normal')
subplot(223)
imagesc(squeeze(abs(mcft_out_ref(end,8,:,:))))
set(gca,'ydir','normal')
subplot(224)
imagesc(squeeze(abs(fftshift(H_ref(end,8,:,:)))))
set(gca,'ydir','normal')

%% inv_mcft

% original function
x_org = inv_mcft(mcft_out_org,cqt_params_out_org,H_org);

% refactored function
x_ref = inv_mcft_refactored(mcft_out_ref,cqt_params_out_ref,H_ref);

% errors
rec_err = norm(x_org - x)/norm(x);
sig_err = norm(x_org - x_ref)/norm(x_org)





