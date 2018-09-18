function [src_est,sdr,sir,sar]=separability_stft(src_trg,mix,stft_params,ibm_thr)

% This function measures the separability of the STFT domain
%
% Inputs:
% src_trg: N*Lt matrix containing N time-domain target sources
% mix: time-domain mixture of N sources
% stft_params: structrue array containing stft parameters including
%              fs,wintype,winlen,ovp
% ibm_thr: 1*m vecotrs containing m ideal binary masking 
%          thresholds (in dB)
%
% Outputs:
% src_est: 1*m cell array containing estimated sources for m ibm thresholds
% sdr: 1*m cell array containing sdr values for sources (each cell constains 
%      N sdr values.
% sir: 1*m cell array containing sir values
% sar: 1*m cell array containing sar values

% dimensions
[num_sources,sig_len] = size(src_trg);

% STFT parameters
fs = stft_params.fs;
wintype = stft_params.wintype;
winlen = stft_params.winlen;
ovp = stft_params.ovp;

% number of threshold values
num_thr = length(ibm_thr);

% compute the STFT of all input signals (mixture & sources)
src_stft = cell(num_sources,1);

for i=1:num_sources
    src_temp = src_trg(i,:);
    [src_stft_temp,~,~,~] = f_stft(src_temp,winlen,wintype,ovp,winlen,fs);
    
    src_stft{i} = src_stft_temp;
end

mix_stft = f_stft(mix,winlen,wintype,ovp,winlen,fs);

% compute source to interference ratios for all sources
src_snr = cell(1,num_sources);
src_num_all = 1:num_sources;

for i=1:num_sources
    % target source
    src_temp = src_stft{i};
    src_temp_mag = abs(src_temp);
    
    % interfering sources
    interf_nums = src_num_all(~ismember(src_num_all,i));
    interf_temp = zeros(size(src_temp));
    for j=1:num_sources-1
        interf_temp = interf_temp + src_stft{interf_nums(j)};
    end
    interf_temp_mag = abs(interf_temp);
        
    % source to interference
    src_snr{i} = 20*log10((src_temp_mag+eps)./(interf_temp_mag+eps));
    
end

% compute ideal binary masked version of all sources
src_masked = cell(num_thr,num_sources);
for i=1:num_thr
    for j=1:num_sources
        
        ibm_temp = src_snr{j}>ibm_thr(i);
        src_masked{i,j} = ibm_temp.*mix_stft;
    end
end

% compute the ISTFT of all estimated sources
src_est = cell(1,num_thr);

for i=1:num_thr
    src_est_temp = zeros(num_sources,sig_len);
    for j=1:num_sources
        est_stft_temp = src_masked{i,j};
        est_time_temp = f_istft(est_stft_temp,winlen,wintype,ovp);
        src_est_temp(j,:) = est_time_temp(1:sig_len);
    end
    src_est{i} = real(src_est_temp);
end

% compute the bss-eval measures
sdr = cell(1,num_thr);
sir = cell(1,num_thr);
sar = cell(1,num_thr);

for i=1:num_thr    
    src_est_temp = src_est{i};
    [sdr{i},sir{i},sar{i},~]=bss_eval_sources(src_est_temp,src_trg);
end


end


