function est_signal = inv_mcft_light(mcft_in,inv_bundle)

% This function reconstructs a time-domain signal given the subsampled
% version of its Multi-resolution Common Fate Transform (MCFT). 
% The intermediary time-frequency domain representation is the 
% Constant-Q Transform (CQT) of the audio signal, which is computed 
% using the invertible and optimized CQT implementation proposed by 
% Schorkhuber et al.:
%
% Toolbox webpage: 
% http://www.cs.tut.fi/sgn/arg/CQT/
% Reference: 
% Schörkhuber et al. "A Matlab toolbox for efficient perfect 
% reconstruction time-frequency transforms with log-frequency resolution."
% Audio Engineering Society Conference: 53rd International Conference: 
% Semantic Audio. Audio Engineering Society, 2014. 
%
% Inputs: 
% mcft_in: cell array, matrix, or structure array containing mcft
%           coefficients
% inv_bundle: structure array containing the information required for the
%             reconstruction of the CQT and then the time-domain signal.
% 
%       cqt_params: structure array containing all CQT properties 
%       fbank: cell array containing the critically sampled scale-rate filter bank
%       filt_ctr_posit: n_scale * n_rate * 2 matrix containig filter centers (low,band,high)
%            in sample # in the scale-rate domain (1st element: scale, 2nd: rate) 
%       filt_range: n_scale * n_rate cell array containing the support of the
%             filter along scale and rate axes in the original sampling
%             grid (transform domain sampled at the actual rate)
%       mcft_scale_len (if output_format = 'matrix'): to be used in mat2cell()
%       mcft_rate_len (if output_format = 'matrix'): to be used in mat2cell()
%
% Ouput:
% est_signal: vector containing the reconstructed time-domain signal
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% Determine the data format of the input MCFT 
 % If the format is 'struct' or 'matrix' convert it to cell array

% total_fbank_size = numel(cell2mat(fbank));

if isstruct(mcft_in) % only zero-padded
        
    mcft_in = struct_to_cell(mcft_in);
    
elseif ~iscell(mcft_in) % is a matrix
    
    mcft_scale_len = inv_bundle.mcft_scale_len;
    mcft_rate_len = inv_bundle.mcft_rate_len;
    
    mcft_in = mat2cell(mcft_in, mcft_scale_len, mcft_rate_len);
        
end

% total_mcft_size = numel(cell2mat(mcft_in));
% 
% if total_mcft_size > total_fbank_size
%     zpadding = 'on';
% else 
%     zpadding = 'off';
% end


%% MCFT to CQT

fbank = inv_bundle.fbank;
filt_ctr_posit = inv_bundle.filt_ctr_posit;
filt_range = inv_bundle.filt_range;
sig_cqt_size = inv_bundle.sig_cqt_size;

disp('Reconstructing the CQT...');
est_sig_cqt = mcft_to_cqt_light(mcft_in,fbank,filt_ctr_posit,filt_range,sig_cqt_size);

%% CQT to time-domain signal

sig_cq_struct = inv_bundle.cqt_params;
sig_cq_struct.c = est_sig_cqt;

disp('Reconstructing the time-domain signal...');
est_signal = icqt(sig_cq_struct);

end


function mcft_cell = struct_to_cell(mcft_struct)

% This function reorders the same-size sections of the filtered CQT, stored
% as 4D matrices in a structure array, and stores them in a cell array
% format.

% Note: only zero-padded version can be stored with a structure format
% Note: lowpass filters are zero-padded to the bandpass size 

% Extract the 4D matrices containing same-size slices
mcft_sband_rband = mcft_struct.sband_rband;
mcft_sband_rhigh = mcft_struct.sband_rhigh;
mcft_shigh_rband = mcft_struct.shigh_rband;
mcft_shigh_rhigh = mcft_struct.shigh_rhigh;

% Convert to cell array slices
[n_sband,n_rband,~,~] = size(mcft_sband_rband);

mcft_sband_rband_cell = cell(n_sband,n_rband);
for i = 1:n_sband
    for j = 1:n_rband
        mcft_sband_rband_cell{i,j} = squeeze(mcft_sband_rband(i,j,:,:));
    end
end

mcft_sband_rhigh_cell = cell(n_sband,2);
for i = 1:n_sband
    mcft_sband_rhigh_cell{i,1} = squeeze(mcft_sband_rhigh(i,1,:,:));
    mcft_sband_rhigh_cell{i,2} = squeeze(mcft_sband_rhigh(i,2,:,:));
end
    
mcft_shigh_rband_cell = cell(2,n_rband);
for i = 1:n_rband
   mcft_shigh_rband_cell{1,i} = squeeze(mcft_shigh_rband(1,i,:,:));
   mcft_shigh_rband_cell{2,i} = squeeze(mcft_shigh_rband(2,i,:,:));
end    

mcft_shigh_rhigh_cell = cell(2,2);
for i = 1:2
    mcft_shigh_rhigh_cell{i,1} = squeeze(mcft_shigh_rhigh(i,1,:,:));
    mcft_shigh_rhigh_cell{i,2} = squeeze(mcft_shigh_rhigh(i,2,:,:));
end


% Concatenate cell array slices
n_scale_ctrs = n_sband + 2;
n_rate_ctrs = n_rband + 2;

scale_high_idx = ceil(n_scale_ctrs/2);
rate_high_idx = ceil(n_rate_ctrs/2);

mcft_cell = cell(n_scale_ctrs, n_rate_ctrs);

mcft_cell(1:scale_high_idx-1, 1:rate_high_idx-1) = ...
    mcft_sband_rband_cell(1:scale_high_idx-1, 1:rate_high_idx-1);
mcft_cell(1:scale_high_idx-1, rate_high_idx+2:end) = ...
    mcft_sband_rband_cell(1:scale_high_idx-1, rate_high_idx:end);
mcft_cell(scale_high_idx+2:end, 1:rate_high_idx-1) = ...
    mcft_sband_rband_cell(scale_high_idx:end, 1:rate_high_idx-1);
mcft_cell(scale_high_idx+2:end, rate_high_idx+2:end) = ...
    mcft_sband_rband_cell(scale_high_idx:end, rate_high_idx:end);

mcft_cell(scale_high_idx:scale_high_idx+1,1:rate_high_idx-1) = ...
    mcft_shigh_rband_cell(:,1:rate_high_idx-1);
mcft_cell(scale_high_idx:scale_high_idx+1,rate_high_idx+2:end) = ...
    mcft_shigh_rband_cell(:,rate_high_idx:end);

mcft_cell(1:scale_high_idx-1,rate_high_idx:rate_high_idx+1) = ...
    mcft_sband_rhigh_cell(1:scale_high_idx-1,:);
mcft_cell(scale_high_idx+2:end,rate_high_idx:rate_high_idx+1) = ...
    mcft_sband_rhigh_cell(scale_high_idx:end,:);

mcft_cell(scale_high_idx:scale_high_idx+1,rate_high_idx:rate_high_idx+1) = ...
    mcft_shigh_rhigh_cell;


end














