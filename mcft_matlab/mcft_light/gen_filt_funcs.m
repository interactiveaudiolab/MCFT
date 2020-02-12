function filt_out = gen_filt_funcs(filt_name,varargin) 

% This function generates the value of filter over a specified range.
% The filter can be a function in the original domain (e.g.,frequency) or
% in the transform domain (e.g., scale)
%
% Inputs: 
% filt_name: character string indicating the name of the filter function
%
% Optional input arguments can be supplied like this:
% 
%    gen_filt_funcs(filt_name,'samp_vec',samp_vec)
%
% The arguments must be character strings followed by the argument value:
%
% 'samp_vec',samp_vec: vector containing the position of samples
% 'filt_len',filt_len: output length (in samples)
% 'dialation',dialation: filter dialation factor 
% 'time_const', time_const: time constant of the exponential term
% 'sigma', sigma: Gabor filter parameter, similar to standard diviation 
% 
% Output:
% filt_out: vector containing the values of the filter
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

% Extract optional input arguments
func_inputs = inputParser;
addParameter(func_inputs,'samp_vec',[],@isnumeric);
addParameter(func_inputs,'filt_len',[],@isnumeric);
addParameter(func_inputs,'dialation',1,@isnumeric);
addParameter(func_inputs,'time_const',1,@isnumeric);
addParameter(func_inputs,'sigma',1,@isnumeric);

parse(func_inputs,varargin{:})
samp_vec = func_inputs.Results.samp_vec;
filt_len = func_inputs.Results.filt_len;
dialation = func_inputs.Results.dialation;
time_const = func_inputs.Results.time_const;
sigma = func_inputs.Results.sigma;

% Check inputs
if nargin < 2
    error('Not enough input arguments!');
elseif isempty(samp_vec) && isempty(filt_len)
    error(['Not enough input arguments: either the vector of ',...
        'sample positions or the filter length must be provided!'])
end

% Generate the sample positions if not provided
if isempty(samp_vec) && ~isempty(filt_len)
    
    % For even filter length the sampling interval is [0:0.5-1/filt_len, -0.5:-1/filt_len]
    if mod(filt_len,2) == 0
        samp_vec = [0 : 1/filt_len : 0.5-1/filt_len, ...
            -0.5 : 1/filt_len : -1/filt_len]';
    
    % For even filter length the sampling interval is [0:0.5-1/(2*filt_len), -0.5+1/(2*filt_len):-2/filt_len]
    else
        samp_vec = [0 : 1/filt_len : 0.5-0.5/filt_len, ...
            -0.5+0.5/filt_len : 1/filt_len : -2/filt_len]';
    end
end

filt_function = gen_func_handle(filt_name);

% Generate filter handles 
switch filt_name    
    case 'hann'
        filt_out = filt_function(samp_vec);
        
    case 'gaborlike' % in the original domain (frequency)
        filt_out = filt_function(samp_vec,dialation);
        
        % make it even so the transform is real
        filt_len = length(filt_out);
        filt_out = [filt_out(1:floor(filt_len/2)+1), filt_out(ceil(filt_len/2):-1:2)];
            
    case 'gaborlike_fourier'
        filt_out = filt_function(samp_vec,dialation);
        
    case 'gabor'
        filt_out = filt_function(samp_vec,dialation,sigma);
        
        % make it even so the transform is real
        filt_len = length(samp_vec);
        filt_out = [filt_out(1:floor(filt_len/2)+1), filt_out(ceil(filt_len/2):-1:2)];
        
        
        % remove dc ???????? 
        
        
    case 'gabor_fourier'
        filt_out = filt_function(samp_vec,dialation,sigma);
        
                % remove dc ???????? 

       
    case 'gammatone' % in the original domain (time)
        filt_out = filt_function(samp_vec,dialation,time_const);
        
        % set the DC value to zero
        filt_out = filt_out - mean(filt_out); 
        
        % make sure the phase is also set to zero 
        % to avoid any computational error (later)
                      
    case 'gammatone_fourier' % in the transform domain (rate)
        filt_out = filt_function(samp_vec,dialation,time_const);
        
                % remove dc ???????? 

        
    case 'damped_sine' 
        filt_out = filt_function(samp_vec,dialation,time_const);
        
        % set the DC value to zero
        filt_out = filt_out - mean(filt_out); 
        
    case 'damped_sine_fourier'
        filt_out = filt_function(samp_vec,dialation,time_const);
        
                % remove dc ???????? 

   
    otherwise
        error('Unknown window function: %s.',name);
end

% Force the window to 0 outside (-.5,.5)
% filt_out = filt_out .* (abs(samp_vec) < .5);    


end


function handle_out = gen_func_handle(filt_name)

% This function generates a function handle for a given filter type
% 
% Inputs:
% filt_name: character string indicating the name of the filter function
%
% Outpus:
% handle_out: function handle 

% Generate the function handle for the given filter name
switch filt_name    
    case 'hann'
        handle_out = @(samp_vec) 0.5 + 0.5 * cos(2 * pi * samp_vec);
        
    case 'gaborlike' % in the original domain (frequency)
        handle_out = @(samp_vec, dialation) dialation * ...
            (1 - 2*(pi * dialation * samp_vec).^2) .* ...
            exp(-((pi * dialation * samp_vec).^2));

    case 'gaborlike_fourier' % in the transform domain (scale)
        handle_out = @(samp_vec, dialation) (samp_vec / dialation).^2 .* ...
            exp(1 - (samp_vec/dialation).^2);
                
    case 'gabor' % in the original domain (frequency)
        handle_out = @(samp_vec,dialation,sigma) cos(2 * pi * dialation * samp_vec) .* ...
            exp(- (dialation * samp_vec).^2 / (2 * sigma ^ 2));
        
    case 'gabor_fourier' % in the transform domain (scale)
        handle_out = @(samp_vec,dialation,sigma) exp(-2 * pi^2 * sigma^2 * (samp_vec/dialation - 1).^2) + ...
            exp(-2 * pi^2 * sigma^2 * (samp_vec/dialation + 1).^2);
        
    case 'gammatone' % in the original domain (time)
        handle_out = @(samp_vec,dialation,time_const) dialation * ...
            (dialation * samp_vec).^2 .* ...
            exp(-time_const * dialation * samp_vec) .* ...
            sin(2 * pi * dialation * samp_vec);
                    
    case 'gammatone_fourier' % in the transform domain (rate)
        term1 = @(samp_vec,dialation,time_const) (time_const + 1j*2*pi*(samp_vec/dialation - 1)).^3;
        term2 = @(samp_vec,dialation,time_const) (time_const + 1j*2*pi*(samp_vec/dialation + 1)).^3;
        handle_out = @(samp_vec,dialation,time_const) 1j*(term1(samp_vec,dialation,time_const) - ...
            term2(samp_vec,dialation,time_const)) ./ (term1(samp_vec,dialation,time_const) .* ...
            term2(samp_vec,dialation,time_const));
    
    case 'damped_sine' % in the original domain (time)
        handle_out = @(samp_vec,dialation,time_const) dialation * ...
            exp(-time_const * dialation * samp_vec) .* ...
            sin(2 * pi * dialation * samp_vec);
            
    case 'damped_sine_fourier' % in the transform domain (rate)
        handle_out = @(samp_vec,dialation,time_const) (2 * pi) ./ ...
            (time_const^2 - 4 * pi^2 * ((samp_vec/dialation).^2 - 1) + ...
            1j * 4 * pi * time_const * (samp_vec/dialation));
        
    otherwise
        error('Unknown window function: %s.',filt_name);
        
end

end

