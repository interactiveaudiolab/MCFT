function ncut_measure=normcut_rep_blk(mix_rep,src_rep,ibm_thr,sim_width,dist,blk_len)

% This function receives constituent sources of a mixture in a
% particular representation domain (STFT, CQT, 2DFT, CFT, MCFT)
% and computes the clusterability measure (normalized cut loss value) as a
% function of ideal binary masking threshold and the width of the
% similarity kernel in the representation domain. 
% Main assumption: the input is high dimensional so it is divided into
% block to avoid exceeding the matrix size limit.
% 
% Inputs:
% mix_rep: (n1*n2*...*nd) matrix containing the mixture
% src_rep: 1*N cell array containing N d-dimensional source
% representations.
% ibm_thr: 1*m vecotrs containing m ideal binary masking thresholds (in dB)
% sim_width: 1*k vector containing k values the similarity kernel width 
%            (in units of the representation domain, 
%              e.g. # of t-f bins for STFT)
% blk_len: (optional) scalar indicating the vlock length in # of samples
%           Default: 100
% dist: (optional) string indicating the distance measure used to compute 
%       the similarity values (see "pdist" for a complete list of 
%       distance measures). Default: 'euclidean'
%
% Ouputs: 
% ncut_measure: an m*k matrix containing normalized cut loss vlues for 
%               m masking thresholds and k similarity kernel widths
%
% Reference: F. R. Bach and M. I. Jordan, 
%            ?Learning spectral clustering, with appli- cation to speech separation,? 
%             Journal of Machine Learning Research, vol. 7, no. Oct, pp. 1963?2001, 2006.

% check inputs
if nargin<4
    dist = 'euclidean';
    blk_len = 100;
elseif nargin<5
    blk_len = 100;
end

% number of parameters
num_thr = length(ibm_thr);
num_sim_width = length(sim_width);

% number of sources
num_sources = length(src_rep);

% compute signal to interference for all sources
src_snr = cell(1,num_sources);
src_num_all = 1:num_sources;

for i=1:num_sources
    
    % target source
    src_temp = src_rep{1,i};
    src_temp_mag = abs(src_temp);
    
    % interfering sources
    interf_nums = src_num_all(~ismember(src_num_all,i));
    interf_temp = zeros(size(src_temp));
    for j=1:num_sources-1
        interf_temp = interf_temp + src_rep{1,interf_nums(j)};
    end
    interf_temp_mag = abs(interf_temp);
    
    % source to interference
    src_snr{i} = 20*log10((src_temp_mag+eps)./(interf_temp_mag+eps));
    
end

% compute sourse estimates by ideal binary masking
mix_mag = abs(mix_rep);

ibms = cell(num_thr,num_sources);

src_masked = cell(num_thr,num_sources);
src_masked_sum = cell(num_thr,1);  % all groups combined
for i=1:num_thr
    src_sum_temp = zeros(size(src_rep{1}));
    for j=1:num_sources
        
        ibm_temp = src_snr{j}>ibm_thr(i);
        ibms{i,j} = ibm_temp;
        
        src_masked_temp = mix_mag.*ibm_temp;
        src_masked_temp = src_masked_temp/max(src_masked_temp(:));
        src_masked_temp = src_masked_temp.*(src_masked_temp>0.1);
        
        src_masked{i,j} = src_masked_temp;
                
        src_sum_temp = src_sum_temp + src_masked_temp;
    end
    src_masked_sum{i} = src_sum_temp;
end

% compute indicator vectors (binary vectors indicating cluster membership)
ind_vecs = cell(num_thr,num_sources);
for i=1:num_thr
    src_sum_temp = src_masked_sum{i};
    % find all nonzero values
    idx_E = find(src_sum_temp);
    for j=1:num_sources
        src_masked_temp = src_masked{i,j};
        % find nonzero values for source j
        idx_e = find(src_masked_temp);
        % form the indicator vector
        ind_vecs{i,j} = ismember(idx_E,idx_e);
    end
end

% find the graph nodes (nonzero elements in the space)
src_dims = size(src_rep{1});
num_src_dims = length(src_dims);

sub_nz = cell(num_thr,1);
for i=1:num_thr
    src_sum_temp = src_masked_sum{i};
    idx_nz = find(src_sum_temp);
    sub_nz_temp = cell(1,num_src_dims);
    [sub_nz_temp{:}] = ind2sub(src_dims,idx_nz);
    sub_nz{i} = cell2mat(sub_nz_temp);
end

% divide samples in each thr group into blocks of length blk_len
sub_nz_blk = cell(num_thr,1);
for i=1:num_thr
    sub_nz_temp = sub_nz{i};
    total_len = size(sub_nz_temp,1);
    if total_len<blk_len
        error('Decrease the block length!')
    end
    
    num_blk = floor(total_len/blk_len);
    disp(['num_blk: ',num2str(num_blk)]);
    
    sub_nz_blk{i} = cell(1,num_blk+1);
    for j=1:num_blk
        samp_range = (j-1)*blk_len+1:j*blk_len;
        sub_nz_blk{i}{j} = sub_nz_temp(samp_range,:);
    end
    last_blk_range = num_blk*blk_len+1:total_len;
    sub_nz_blk{i}{num_blk+1} = sub_nz_temp(last_blk_range,:);
end

% compute clusterability measures
disp('Computing n-cut ...');

ncut_measure = zeros(num_thr,num_sim_width);

for i=1:num_thr
    
    sub_nz_temp = sub_nz{i};    
    
    total_len = size(sub_nz_temp,1); 
    num_blk = floor(total_len/blk_len);

    ones_col = ones(total_len,1);

    sub_nz_blk_temp = sub_nz_blk{i};
    ind_vecs_temp = cell2mat(ind_vecs(i,:));
        

    for j=1:num_sim_width
        disp(['thr: ',num2str(ibm_thr(i)),' sim_width: ',num2str(sim_width(j))]);

        % process by block                
        D_ind_cell = cell(num_blk+1,1);
        L_ind_cell = cell(num_blk+1,1);   
                
        for k=1:num_blk+1
            disp(['sim: ',num2str(j),' - blk: ',num2str(k),'/',num2str(num_blk)])
            % compute pairwise distances 
            blk_pdists = pdist2(sub_nz_blk_temp{k},sub_nz_temp,dist);
            
            % compute similarity values
            if strcmp(dist,'euclidean')
                blk_sim = exp(-blk_pdists.^2/sim_width(j)^2);
            else
                blk_sim = exp(-blk_pdists/sim_width(j));
            end

            % compute the D and L matrices
            blk_D = zeros(size(blk_sim));
            blk_start = (k-1)*blk_len+1;
            blk_end = min(k*blk_len, total_len);

            blk_D(:,blk_start:blk_end) = diag(blk_sim*ones_col);

            % compute the Laplacian matrix
            blk_L = blk_D - blk_sim;

            % compute D*e and L*e
            D_ind_cell{k} = blk_D * ind_vecs_temp;
            L_ind_cell{k} = blk_L * ind_vecs_temp;
            
        end
        
        % extract D and L values from cell arrays
        D_ind_vec = cell2mat(D_ind_cell);
        L_ind_vec = cell2mat(L_ind_cell);
        
                
        ncut_numerator = diag(ind_vecs_temp' * L_ind_vec)+eps;
        ncut_denominator = diag(ind_vecs_temp' * D_ind_vec)+eps;
        ncut_measure(i,j) = sum(ncut_numerator./ncut_denominator);
        
    end
end



end



