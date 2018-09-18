function loss =clusterability_stft(src_trg,mix,stft_params,ncut_params)

% This function measures the clusterability of the STFT domain
%
% Inputs:
% src_trg: N*Lt matrix containing N time-domain target sources
% mix: time-domain mixture of N sources
% stft_params: structrue array containing stft parameters including
%              fs,wintype,winlen,ovp
% ncut_params: structure array containing normalized-cut parameters 
%              including: ibm_thr,sim_width,dist,blk_len
%
% Output:
% loss: an m*k matrix containing normalized cut loss vlues for 
%               m masking thresholds and k similarity kernel widths

% number of sources
num_sources = size(src_trg,1);

% STFT parameters
fs = stft_params.fs;
wintype = stft_params.wintype;
winlen = stft_params.winlen;
ovp = stft_params.ovp;

% ncut parameters
ibm_thr = ncut_params.ibm_thr;
sim_width = ncut_params.sim_width;
dist = ncut_params.dist;
blk_len = ncut_params.blk_len;

% compute the STFT of all input signals (mixture & sources)
src_stft = cell(1,num_sources);

for i=1:num_sources
    src_temp = src_trg(i,:);
    [src_stft_temp,~,~,~] = f_stft(src_temp,winlen,wintype,ovp,winlen,fs);
    
    src_stft{i} = src_stft_temp;
end

mix_stft = f_stft(mix,winlen,wintype,ovp,winlen,fs);

% compute clusterability measure
loss = normcut_rep_blk(mix_stft,src_stft,ibm_thr,sim_width,dist,blk_len);

end

