%% This script creates a filterbank class 
 % for the Multi-resolution Common Fate Transform 
 
 classdef MCFTfilterbank
     % What is this class? 
     % What are its properties?
     % What are its methods?
     
     %% Properties of the filterbank:
      % 
      %%% General filterbank properties:
      % fbank_domain: character string indicating whether the filter should
      %               be computed in the original (time-frequency:'tf') domain
      %               or in the transform (scale-rate:'sr') domain
      % 
      % complex_specgram: n_freq by n_time matrix containing a 2d 
      %                   complex spectrogram. If provided, the filterbank 
      %                   will be modulated with the phase of the spectrogram. 
      %                   Otherwise, the function the original (unmodulated)
      %                   set of filters will will be returned   
      % 
      %%% Scale filterbank properties:
      % scale_filt_func: character string specifying the name of the function 
      %                    used for computing the values of the scale filters, 
      %                    e.g., 'Gabor', or 'Hann'
      % scale_ctrs: vector containing filter centers along the scale axis
      % spec_samprate: sample rate along the freq. axis (cyc/oct)
      % nfft_scale: number of scale fft points
      % sigma: Gabor filter parameter, similar to standard diviation
      % scale_res: number of bins per octave on the scale axis, default = 1  
      % scale_max: the center of the highest bandpass filter
      %                  default: 2^(nextpow2(samprate_spec/2)-1) 
      %                  (last power2 number before the Nyquist rate)
      %                  e.g., scale_max = 8 if samp rate = 24 
      % scale_min: the center of the lowest bandpass filter
      %                  default: 2^0   
      %
      %%% Rate filterbank properties:
      % rate_filt_func: character string specifying the name of the function 
      %                    used for computing the values of the rate filters, 
      %                    e.g., 'Gabor', or 'Hann'
      % rate_ctrs: vector containing filter centers along the rate axis
      % temp_samprate: sample rate along the time axis
      % nfft_rate: number of rate fft points
      % time_const: exponent coefficient of the exponential term
      % rate_res: number of bins per octave on the rate axis, default = 1
      % rate_max: the center of the highest pandpass filter
      %                  defa ult: 2^(nextpow2(samprate_temp/2)-1) 
      %                  (last power2 number before the Nyquist rate)
      % rate_min: the center of the lowest filter
      %                  default: 2^0
       
     properties
         % General filterbank properties
         fbank_domain = 'sr';            % default: 'sr' (for subsampling)
         complex_specgram = [];          % default: []
         
         % Scale filterbank properties
         scale_filt_func = 'gabor';      % default: 'gabor'
         scale_ctrs = [];                % default: []
         spec_samprate = [];             % default: []
         nfft_scale = [];                % default: []
         sigma = 1;                      % default: 1
         scale_res = 1                   % default: 1
         scale_max = [];                 % default: [] 2^(nextpow2(spec_samprate/2)-1)
         scale_min = 2^0;                % default: 2^0
         
         % Rate filterbank properties
         rate_filt_func = 'gammatone';   % default: 'gammatone'
         rate_ctrs = [];                 % default: []
         temp_samprate = [];             % default: []
         nfft_rate = [];                 % default: []
         time_const = 1;                 % default: 1
         rate_res = 1;                   % default: 1
         rate_max = [];                  % default: [] 2^(nextpow2(temp_samprate/2)-1)
         rate_min = 2^0;                 % default: 2^0
         
     end
     
     
     
     %% Methods of the filterbank
     methods
         
         %%%  Constructor method %%%
         % Construct (initialize) the filterbank object with the initial
         % property values, which can be specified as Name-Value pairs
         % e.g., fbank_obj = MCFTfilterbank('scale_ctrs',scale_ctrs)
         
         function obj = MCFTfilterbank(varargin)
             func_inputs = inputParser;
             
             % Scale filterbank parameters:
             addParameter(func_inputs,'spec_samprate',[],@isnumeric);
             addParameter(func_inputs,'nfft_scale',[],@isnumeric);
             addParameter(func_inputs,'scale_res',1,@isnumeric);
             
             % Rate filterbank parameters
             addParameter(func_inputs,'temp_samprate',[],@isnumeric);
             addParameter(func_inputs,'nfft_rate',[],@isnumeric);
             addParameter(func_inputs,'time_const',1,@isnumeric);
             addParameter(func_inputs,'rate_res',1,@isnumeric);
             
             % Parse the inputs
             parse(func_inputs,varargin{:})    
                          
             obj.spec_samprate = func_inputs.Results.spec_samprate;
             obj.nfft_scale = func_inputs.Results.nfft_scale;
             obj.scale_res = func_inputs.Results.scale_res;
             
             obj.temp_samprate = func_inputs.Results.temp_samprate;
             obj.nfft_rate = func_inputs.Results.nfft_rate;
             obj.time_const = func_inputs.Results.time_const;
             obj.rate_res = func_inputs.Results.rate_res;
             
             
             % compute initial scale_max value
             if isempty(obj.scale_max) && ~isempty(obj.spec_samprate)
                 obj.scale_max = 2^(nextpow2(obj.spec_samprate/2)-1);
             end
             
             % compute initial scale_max value
             if isempty(obj.rate_max) && ~isempty(obj.temp_samprate)
                 obj.rate_max = 2^(nextpow2(obj.temp_samprate/2)-1);
             end
             
         end
         
         %%% Generate filter centers
         % Scale centers
         function scale_ctrs = gen_scale_ctrs(obj)
               scale_params = struct('filt_type','scale','filt_res',obj.scale_res,...
                   'filt_nfft',obj.nfft_scale,'samprate',obj.spec_samprate);
               scale_ctrs = filt_default_centers_light(scale_params);
         end
         
         % Rate centers      
         function rate_ctrs = gen_rate_ctrs(obj)
               rate_params = struct('filt_type','rate','filt_res',obj.rate_res,...
                   'filt_nfft',obj.nfft_rate,'samprate',obj.temp_samprate);
               rate_ctrs = filt_default_centers_light(rate_params);
         end
         
         
         %%% Filter function generator
         function filt_out = gen_filt_scale(obj,scale_ctr)
             scale_filt_len = obj.nfft_scale;
             scale_vec = obj.spec_samprate * (-0.5:1/scale_filt_len:0.5-1/scale_filt_len); 
             
             filt_out = gen_filt_funcs('gabor_fourier','samp_vec',scale_vec,...
                 'dialation',scale_ctr,'sigma',obj.sigma);       
         end
         
         
        
         
     end
     
 end