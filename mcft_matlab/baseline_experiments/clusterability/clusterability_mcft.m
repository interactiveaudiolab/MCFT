function loss =clusterability_mcft(src_trg,mix,cqt_params,ncut_params,beta)

% This function measures the clusterability of the MCFT domain
%
% Inputs:
% src_trg: N*Lt matrix containing N time-domain target sources
% mix: time-domain mixture of N sources
% cqt_params: structrue array containing stft parameters including
%              fs,fmin,fmax,fres,gamma
% ncut_params: structure array containing normalized-cut parameters 
%              including: ibm_thr,sim_width,dist,blk_len
% beta: (optional) filterbank time constant. default: 1
%
% Output:
% loss: an m*k matrix containing normalized cut loss vlues for 
%               m masking thresholds and k similarity kernel widths

if nargin<5
    beta = 1;
end

% dimensions
[num_sources,sig_len] = size(src_trg);

% CQT parameters
fs = cqt_params.fs;
fmin = cqt_params.fmin;
fmax = cqt_params.fmax;
fres = cqt_params.fres;
gamma = cqt_params.gamma;

% ncut parameters
ibm_thr = ncut_params.ibm_thr;
sim_width = ncut_params.sim_width;
dist = ncut_params.dist;
blk_len = ncut_params.blk_len;

% compute the CQT of all input signals (mixture & sources)
src_cqt_struct = cell(1,num_sources);

for i=1:num_sources
    src_temp = src_trg(i,:);
    src_cqt_temp = cqt(src_temp,fres,fs,fmin,fmax,'rasterize','full','gamma',gamma);
    src_cqt_struct{i} = src_cqt_temp;
end

mix_cqt_struct = cqt(mix,fres,fs,fmin,fmax,'rasterize','full','gamma',gamma);
mix_cqt = mix_cqt_struct.c;

% compute the filter bank
[num_f,num_t] = size(mix_cqt);
nfft_s = num_f;
nfft_r = num_t;
samprate_s = fres;
sig_dur = sig_len/fs;
samprate_r = floor(nfft_r/sig_dur);

[sv,rv]=filt_default_centers(nfft_s,nfft_r,samprate_s,samprate_r);

h_params=struct('samprate_spec',samprate_s,'samprate_temp',samprate_r,'time_const',beta);
[~,H]=gen_fbank_hsr(sv,rv,nfft_s,nfft_r,h_params,mix_cqt); 

% compute the MCFT of all input signals (mixture & sources)
src_mcft = cell(1,num_sources);

for i=1:num_sources
    src_cqt_temp = src_cqt_struct{i}.c;
    src_mcft_temp = cqt_to_mcft(src_cqt_temp,H);
    src_mcft{i} = src_mcft_temp;
end

mix_mcft = cqt_to_mcft(mix_cqt,H);

% compute clusterability measure
loss = normcut_rep_blk(mix_mcft,src_mcft,ibm_thr,sim_width,dist,blk_len);


end