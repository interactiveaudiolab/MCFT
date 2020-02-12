function [fbank,filt_ctrs,ctr_shift,filt_len] = gen_rate_fbank_light(ctr_min,ctr_max,...
    filt_res,samprate,nfft,varargin)

% This function generates a multirate scale/rate filterbank.
% Bandpass filters are critically sampled.
%
% Inputs:
% ctr_min: center of the lowest bandpass rate filter
% ctr_max: center of the highest bandpass rate filter
% filt_res: number of bins (filters) per octave
% samprate: temporal sampling rate
% nfft: number of samples in the transform domain
%
% Outputs: 
% fbank: cell array containing constant-Q/variabale-Q rate filters
% ctr_shift: vector containing shifts between the center frequencies
% filt_len: vector containing filter lengths
%
% Optional input arguments can be provided like this:
% 
%   gen_rate_fbank_light(ctr_min,ctr_max,filt_res,samprate,nfft,...
%                      'min_filt_len',min_filt_len)
%
% The optional arguments much be names(character strings) followed by values:
% 
% 'filt_name': character string specifying the name of the filter function
% 'min_filt_len': minimum filter length 
% 'bw_offset': bandwidth offset. If bw_offset = 0 the filterbank is constant-Q.
%        bandwidth = 1/Q_fator * ctr + bw_offset (Q is determined by #bins per
%        octave)
%        bw_offset > 0 improves the resolution of lower filters 
% 'time_const': time constant of the exponential term
%

%%% Extract optional input arguments
func_inputs = inputParser;
addParameter(func_inputs,'filt_name','gammatone_fourier',@ischar);
addParameter(func_inputs,'min_filt_len',4,@isnumeric);
addParameter(func_inputs,'bw_offset',0,@isnumeric);
addParameter(func_inputs,'time_const',1,@isnumeric);

parse(func_inputs,varargin{:})
filt_name = func_inputs.Results.filt_name;
min_filt_len = func_inputs.Results.min_filt_len;
bw_offset = func_inputs.Results.bw_offset;
time_const = func_inputs.Results.time_const;

%%% Compute filter centers and bandwidths
nyq_rate = samprate / 2; % Nyquiest rate

if ctr_max > nyq_rate, ctr_max = nyq_rate; end

fft_res = samprate / nfft;

% index of the highest filter center in log2 scale (strating at 2^0)
ctr_max_idx = floor(filt_res * log2(ctr_max / ctr_min));

% center frequencies of bandpass filters
filt_ctrs = ctr_min .* 2.^((0:ctr_max_idx).' ./ filt_res);

Q_factor = 2^(1/filt_res) - 2^(-1/filt_res);

% bandwidths of bandpass filters
filt_bws = Q_factor * filt_ctrs + bw_offset;

% make sure the support of the highest filter will not exceed nyq_rate
temp_idx = find(filt_ctrs + filt_bws/2 > nyq_rate, 1, 'first');
if ~isempty(temp_idx)
    filt_ctrs = filt_ctrs(1:temp_idx - 1);
    filt_bws = filt_bws(1:temp_idx - 1);
end

% make sure the support of the lowest filter will not exceed DC
temp_idx = find(filt_ctrs - filt_bws/2 < 0, 1, 'last');
if ~isempty(temp_idx)
    filt_ctrs = filt_ctrs(temp_idx + 1 : end);
    filt_bws = filt_bws(temp_idx + 1 : end);
    warning(['ctr_min set to ',num2str(fft_res * floor(filt_ctrs(1)/fft_res),6),'!']);
end
    
% total number of bandpass filters
n_bpass = length(filt_ctrs);

% include DC
% filt_ctrs = [0; filt_ctrs];
filt_ctrs = [fft_res/2; filt_ctrs];

% iclude Nyquist
filt_ctrs = [filt_ctrs; nyq_rate];

% include mirrored centers and bandwidths in (pi,2pi)
filt_ctrs(n_bpass+3 : 2*(n_bpass+1)) = samprate - filt_ctrs(n_bpass+1 : -1 : 2);
filt_bws = [2*ctr_min ; filt_bws ; filt_ctrs(n_bpass+3)-filt_ctrs(n_bpass+1); filt_bws(end:-1:1)];

% filter centers and bandwidths in terms of bin number
filt_ctrs_samp = filt_ctrs / fft_res;
filt_bws_samp = round(filt_bws / fft_res);

% position of filter centers in terms of bin number
ctr_posit = zeros(size(filt_ctrs_samp));

% floor of positive filter centers (round toward zero)
ctr_posit(1:n_bpass+2) = floor(filt_ctrs_samp(1:n_bpass+2));

% ceiling of negative filter centers (round toward zero)
ctr_posit(n_bpass+3:end) = ceil(filt_ctrs_samp(n_bpass+3:end));

% shift in filter centers
% mod(-ctr_posit(end),nfft) is the distance between the last center freq and 2pi (or 0 to first center freq)
ctr_shift = [mod(-ctr_posit(end),nfft) ; diff(ctr_posit)];
    
% map the bandwithds and filter lengths to the next even number
% filt_bws_samp = filt_bws_samp + mod(filt_bt_samp, 2);

% check if filter lengths are shorter than allowed
filt_bws_samp(filt_bws_samp < min_filt_len) = min_filt_len;
filt_len = filt_bws_samp;

% make an array of filter functions
fbank = cell(2*(n_bpass+1), 1);
for i = 1:2*(n_bpass+1)
    len_temp = filt_len(i);
    
     if mod(len_temp,2) == 0 % For even N the sampling interval is [-.5,.5-1/N]
        % r_ctr = 1, 33dB bandwidth -> Q = 1.09
        samp_vec_temp = (2^1 - 2^0) * [0:1/len_temp:.5-1/len_temp,...
                                  -.5:1/len_temp:-1/len_temp]';
                                   
    else % For odd N the sampling interval is [-.5+1/(2N),.5-1/(2N)]
        samp_vec_temp = (2^1 - 2^0) * [0:1/len_temp:.5-.5/len_temp,...
                                  -.5+.5/len_temp:1/len_temp:-1/len_temp]';
    end
    
    % shift and scale the s-axis
    samp_vec_shift = (samp_vec_temp(:) + 1);
    
    % compute the filter 
    filt_temp = gen_filt_funcs(filt_name,'samp_vec',samp_vec_shift,...
        'dialation',1,'time_const',time_const);
    
    if i < n_bpass+3 
        fbank{i} = filt_temp;
    else
        fbank{i} = conj(filt_temp(end:-1:1));
    end
end


% Setup Tukey window for 0- and Nyquist-frequency
for i = [1,n_bpass+2]
    len_temp = filt_len(i);
    len_next = filt_len(i+1);
    
    if len_temp > len_next
        
        filt_mag_temp = ones(len_temp,1);
        samp_vec_temp = [0:1/len_next:.5-1/len_next,-.5:1/len_next:-1/len_next]';
        
        hann_temp = .5 + .5 * cos(2 * pi * samp_vec_temp);
        
        filt_mag_temp((floor(len_temp/2)-floor(len_next/2)+1):(floor(len_temp/2) + ...
            ceil(len_next/2))) = hann_temp;
        
        fbank{i} = filt_mag_temp .* exp(1j*angle(fbank{i}));
%         fbank{i} = fbank{i} / sqrt(len_temp);
    end
end


end








