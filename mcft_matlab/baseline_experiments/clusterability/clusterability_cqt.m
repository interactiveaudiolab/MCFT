function loss =clusterability_cqt(src_trg,mix,cqt_params,ncut_params)

% This function measures the clusterability of the CQT domain
%
% Inputs:
% src_trg: N*Lt matrix containing N time-domain target sources
% mix: time-domain mixture of N sources
% cqt_params: structrue array containing stft parameters including
%              fs,fmin,fmax,fres,gamma
% ncut_params: structure array containing normalized-cut parameters 
%              including: ibm_thr,sim_width,dist,blk_len
%
% Output:
% loss: an m*k matrix containing normalized cut loss vlues for 
%               m masking thresholds and k similarity kernel widths

% number of sources
num_sources = size(src_trg,1);

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
src_cqt = cell(1,num_sources);

for i=1:num_sources
    src_temp = src_trg(i,:);
    src_cqt_temp = cqt(src_temp,fres,fs,fmin,fmax,'rasterize','full','gamma',gamma);
    src_cqt{i} = src_cqt_temp.c;
end

mix_cqt_struct = cqt(mix,fres,fs,fmin,fmax,'rasterize','full','gamma',gamma);
mix_cqt = mix_cqt_struct.c;

% compute clusterability measure
loss = normcut_rep_blk(mix_cqt,src_cqt,ibm_thr,sim_width,dist,blk_len);

end





