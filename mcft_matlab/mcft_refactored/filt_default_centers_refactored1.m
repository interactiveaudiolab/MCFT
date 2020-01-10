function [scale_ctrs,rate_ctrs] = filt_default_centers_refactored1(scale_params,rate_params)

% This function computes the default set of filter centers along scale
% and rate axes. Two cases are considered:
% 1. Inputs only include the resolution, number of fft points, and
%    sample rate. In this case, default values will be used for all
%    filters.
% 2. Inputs also include the minimum and maximum values for bandpass filter
%    centers. In this case, default values will be used for lowpass and 
%    highpass filters.
%
% Inputs:
% scale_params: structure array containing parameters of scale filters:
%               scale_res: number of bins per octvave on the scale axis
%               scale_nfft: number of fft points on the scale axis
%               samprate_spec: sampling rate of the spectral filter
%                             (in cycles per octave)
%               scale_max (optional): the center of the highest bandpass filter
%                          default: scale_max = samprate_spec/2
%               scale_min (optional): the center of the lowest bandpass filter
%                          default: 2^0
%
% rate_params: structure array containing parameters of rate filters:
%              rate_res: number of bins per octave on the rate axis
%              rate_nfft: number of fft points on theh rate axis
%              samprate_temp: sampling rate of the temporal filters 
%                            (in cycles per second)
%              rate_max (optional): the center of the highest bandpass filter
%                                   default: samprate_rate/2
%              rate_min (optional): the center of the lowest bandpass filter
%                        default: 2^0
%
% Outputs: 
% scale_ctrs: vector containing filter centers (along scale axis)
% rate_ctrs: vector containing filter centers (along rate axis)
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

% extract scale parameters
scale_res = scale_params.scale_res;
scale_nfft = scale_params.scale_nfft;
samprate_spec = scale_params.samprate_spec;
scale_max = [];
scale_min = [];
if isfield(scale_params,'scale_max'), scale_max = scale_params.scale_max; end
if isfield(scale_params,'scale_min'), scale_min = scale_params.scale_min; end

% compute scale filter centers
scale_ctrs = filt_centers('scale',scale_res,scale_nfft,samprate_spec,...
    'ctr_max',scale_max,'ctr_min',scale_min);

% extract rate parameters
rate_res = rate_params.rate_res;
rate_nfft = rate_params.rate_nfft;
samprate_temp = rate_params.samprate_temp;
rate_max = [];
rate_min = [];
if isfield(rate_params,'rate_max'), rate_max = rate_params.rate_max; end
if isfield(rate_params,'rate_min'), rate_min = rate_params.rate_min; end

% compute rate filter centers
rate_ctrs = filt_centers('rate',rate_res,rate_nfft,samprate_temp,...
'ctr_max',rate_max,'ctr_min',rate_min);

end


function filt_ctrs = filt_centers(filt_type,bins_per_oct,nfft,samprate,varargin)

% This function computes scale/rate filter centers given transform-domain
% parameters.
%
% Inputs:
% filt_type: string, 'scale' or 'rate' (highpass filter computed
%           differently based on filter type)
% bins_per_oct: number of scale/rate filters per octave
% nfft: number of frequencies of analysis in the scale/rate domain
% samprate: sampling rate in the spectral/temporal domain
% 
% Optional input arguments can be supplied like this:
%
%    filt_centers(filt_type,bins_per_oct,nfft,'ctr_max',ctr_max)
%
% The arguments must be character strings followed by the argument value:
%
% 'ctr_max',ctr_max: the center of the highest bandpass filter
%                    default: 2^(nextpow2(samprate/2)-1) 
%                            (last power2 number before the Nyquist rate)
%                            e.g., ctr_max = 8 if samprate = 24 
% 'ctr_min',ctr_min: the center of the lowest filter
%                            default: 2^0
%
% Outputs:
% filt_ctrs: numpy array containing filter centers

% extract optional inputs
func_inputs = inputParser;
addParameter(func_inputs,'ctr_min',[],@isnumeric);
addParameter(func_inputs,'ctr_max',[],@isnumeric);
parse(func_inputs,varargin{:})

ctr_max = func_inputs.Results.ctr_max;
ctr_min = func_inputs.Results.ctr_min;

if ~isempty(ctr_min) 
    log2_ctr_band_min = log2(ctr_min);    
else
    log2_ctr_band_min = 0;
end

if ~isempty(ctr_max)
    log2_ctr_band_max = log2(ctr_max);
else
    % take the largest center that is smaller than the Nyquist rate
%     log2_ctr_band_max = max(log2_ctr_band_min:1/bins_per_oct:log2(samprate/2));
      log2_ctr_band_max = nextpow2(samprate/2) - 1;
end

% spacing between the frequencies of analysis
grid_res = samprate/nfft;
   
% lowpass filter center
ctr_low = grid_res/2; 

% add 1 to power of 2 if the smallest bandpass center is 
% smaller than the lowpass center
if 2 ^ log2_ctr_band_min <= ctr_low
    log2_ctr_band_min = log2(ctr_low) + 1 ;
end

ctr_band = 2.^(log2_ctr_band_min:1/bins_per_oct:log2_ctr_band_max);

% highpass filter centers
if strcmp(filt_type, 'scale')
    ctr_high = (ctr_band(end) + 3*samprate / 2) / 4;
elseif strcmp(filt_type, 'rate')
    ctr_high = (ctr_band(end) + 1*samprate / 2) / 2;
end

% shift filter centers to the nearest frequency of analysis
ctr_band = round(ctr_band/grid_res) * grid_res;
ctr_high = round(ctr_high/grid_res) * grid_res;

% concatenate all centeres into one vector
ctr_all = [ctr_low, ctr_band, ctr_high];

% remove repeated values (sometimes happens due to rounding
filt_ctrs = unique(ctr_all);

end

