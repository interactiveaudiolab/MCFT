function [filt_tf_domain,filt_sr_domain] = gen_filt_scale_rate_light(scale_ctr,rate_ctr,scale_params,rate_params,filt_dir)

% This function generates a 2D-impulse response in the time-frequncy
% domain with dilation factors S and R: h(omega,tau;S,R)
%
% Inputs:
% scale_ctr: filter center along the scale axis
% rate_ctr: filter center along the rate axis  
% scale_params:structure array containing the parameters of the spectral filter
%          spec_filt_len: length of the spectral filter impulse response
%          spec_samprate: sample rate along the freq. axis (cyc/oct)
%          spec_filt_type: character string argument indicating the filter type
%               ('bandpass','lowpass','highpass')
% rate_params:structure array containing the parameters of the temporal filter
%          time_const: exponent coefficient of the exponential term
%          temp_filt_len: length of the temporal filter impulse response
%          temp_samprate: sample rate along the time axis
%          temp_filt_type: string argument indicating the type of the filter
%                ('bandpass','lowpass','highpass')
% filt_dir: character string type, determines the moving direction of the filter 
%            'none' (full s-r domain)
%            'up' (upward analytic, nonzero over upper left (2nd) and lower right (4th) quadrants)
%            'down' (downward analytic, nonzero over upper right (1st) and lower left (3rd) quadrants)
%
% Outputs
% filt_tf_domain: impulse response of the 2D filter (in time-frequency domain)
% filt_sr_domain: 2D-Fourier transform of the filter impulse response 
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% extract filter parameters

% scale/spectral filter
spec_filt_len = scale_params.spec_filt_len; % length of the spectral filter
spec_samprate = scale_params.spec_samprate; % spectral ripple frequency
spec_filt_type = scale_params.spec_filt_type; % type of the spectral filter

% rate/temporal filter 
beta = rate_params.time_const; % time constant of the temporal filter
temp_filt_len = rate_params.temp_filt_len;  % length of the temporal filter
temp_samprate = rate_params.temp_samprate; % # of frames per second
temp_filt_type = rate_params.temp_filt_type;  % type of the temporal filter


%% frequency and time vectors

% zero-pad filters to the next even number 
spec_filt_len = spec_filt_len + mod(spec_filt_len,2);
temp_filt_len = temp_filt_len + mod(temp_filt_len,2);
    
% frequency and time vectors
freq_vec = (0:spec_filt_len-1)'/spec_samprate;
time_vec = (0:temp_filt_len-1)'/temp_samprate;

%% impulse response of the original scale filter: Gaussian 

spec_filt = scale_ctr * (1 - 2 * (scale_ctr * pi * freq_vec).^2) .* ...
    exp(-((scale_ctr * freq_vec * pi).^2));

% make it even so the transform is real
spec_filt = [spec_filt(1:spec_filt_len/2+1) ; spec_filt(spec_filt_len/2:-1:2)]; 

%% impulse response of the original rate filter    

temp_filt = rate_ctr * (rate_ctr * time_vec).^2 .* exp(-time_vec * beta * rate_ctr)...
    .* sin(2 * pi * rate_ctr * time_vec);
temp_filt = temp_filt - mean(temp_filt); 

% if the magnitude of dc element is set to zero by
% subtracting the mean of hr, make sure the phase is
% also set to zero to avoid any computational error
if abs(mean(temp_filt))<eps
    correct_rate_filt_phase = 1;
else
    correct_rate_filt_phase = 0;
end

%% scale response (Fourier transform of the scale impulse response)

% band-pass scale filter
scale_filt = abs(fft(spec_filt, spec_filt_len)); % discard the negligible imaginary part

% low/high-pass scale filter
if ~strcmp(spec_filt_type, 'bandpass')    
    scale_filt_1 = scale_filt(1:spec_filt_len/2+1, :); 
    scale_filt_1 = scale_filt_1/max(scale_filt_1);
    
    scale_filt_2 = scale_filt(spec_filt_len/2+2:end, :); 
    scale_filt_2 = scale_filt_2/max(scale_filt_2);
    
    [~, max_idx1] = max(scale_filt_1);
    [~, max_idx2] = max(scale_filt_2);
    
    if strcmp(spec_filt_type, 'lowpass')
        scale_filt_1(1:max_idx1-1) = 1;  
        scale_filt_2(max_idx2+1:end) = 1;
        
    elseif strcmp(spec_filt_type,'highpass')
        scale_filt_1(max_idx1+1:end) = 1;
        scale_filt_2(1:max_idx2-1) = 1;
        
    end
    
    scale_filt = [scale_filt_1 ; scale_filt_2]; % form the full magnitude spectrum
end 


%% rate response (Fourier transform of the rate impulse response)

% band-pass rate filter
rate_filt = fft(temp_filt, temp_filt_len); % rate response is complex

% low/high-pass rate filter
if ~strcmp(temp_filt_type, 'bandpass')
    rate_filt_ph = angle(rate_filt); 
    if correct_rate_filt_phase
     rate_filt_ph(1) = 0; % dc element is zero so the phase shouldn't matter
    end
    rate_filt_mag = abs(rate_filt);
    
    rate_filt_mag_1 = rate_filt_mag(1:temp_filt_len/2+1, :); 
    rate_filt_mag_1 = rate_filt_mag_1/max(rate_filt_mag_1);
    
    rate_filt_mag_2 = rate_filt_mag(temp_filt_len/2+2:end, :); 
    rate_filt_mag_2 = rate_filt_mag_2/max(rate_filt_mag_2);
    
    [~, max_idx1] = max(rate_filt_mag_1);
    [~, max_idx2] = max(rate_filt_mag_2);
    
    if strcmp(temp_filt_type, 'lowpass')
       rate_filt_mag_1(1:max_idx1-1) = 1;    
       rate_filt_mag_2(max_idx2+1:end) = 1;
       
    elseif strcmp(temp_filt_type, 'highpass')
       rate_filt_mag_1(max_idx1+1:end) = 1;
       rate_filt_mag_2(1:max_idx2-1) = 1;
       
    end
    
    rate_filt_mag = [rate_filt_mag_1 ; rate_filt_mag_2]; % form the full magnitude spectrum
    rate_filt = rate_filt_mag .* exp(1j * rate_filt_ph); % form the full Fourier transform
end    
 
%% full scale-rate impulse and transform responses
 % filt_sr_full is quadrant separable

filt_sr_full = scale_filt * rate_filt.';

% normalize the filter magnitude:
filt_sr_full_mag = abs(filt_sr_full);
filt_sr_full_mag = filt_sr_full_mag / max(filt_sr_full_mag(:));
filt_sr_full = filt_sr_full_mag .* exp(1j*angle(filt_sr_full));

if strcmp(filt_dir,'up')
    
  % compute the upward version of the scale-rate response
  filt_sr_up = filt_sr_full;
  filt_sr_up(2:spec_filt_len/2+1 , 2:temp_filt_len/2+1) = 0; 
  %filt_sr_up(1:spec_filt_len/2 , 1:temp_filt_len/2)=0; 
  filt_sr_up(spec_filt_len/2+2:end , temp_filt_len/2+2:end) = 0; 
  filt_sr_domain = filt_sr_up;
  
elseif strcmp(filt_dir,'down') % downward ripple
    
  filt_sr_down = filt_sr_full;
  %filt_sr_down(1:spec_filt_len/2 , temp_filt_len/2+2:end)=0;
  %filt_sr_down(spec_filt_len/2+2:end , 1:temp_filt_len/2)=0;
  filt_sr_down(2:spec_filt_len/2+1 , temp_filt_len/2+2:end) = 0;
  filt_sr_down(spec_filt_len/2+2:end , 2:temp_filt_len/2+1) = 0;
  filt_sr_domain = filt_sr_down;
  
else   
   filt_sr_domain = filt_sr_full;
end

filt_tf_domain = ifft2(filt_sr_domain);

if max(imag(filt_tf_domain(:)))<1e-8
    filt_tf_domain = real(filt_tf_domain); 
end

end


