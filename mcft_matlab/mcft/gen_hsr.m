function [hsr,Hsr] = gen_hsr(scale_ctr,rate_ctr,scale_params,rate_params,filt_dir)

% This function generates a 2D-impulse response in the time-frequncy
% domain with dilation factors S and R: h(omega,tau;S,R)
%
% Inputs:
% scale_ctr: filter center along the scale axis
% rate_ctr: filter center along the rate axis  
% scale_params:structure array containing the parameters of the spectral filter
%          hslen: length of the spectral filter impulse response
%          samprate_spec: sample rate along the freq. axis (cyc/oct)
%          type: string argument indicating the filter type
%               ('bandpass','lowpass','highpass')
% rate_params:structure array containing the parameters of the temporal filter
%          time_const: exponent coefficient of the exponential term
%          hrlen: length of the temporal filter impulse response
%          samprate_temp: sample rate along the time axis
%          type: string argument indicating the type of the filter
%                ('bandpass','lowpass','highpass')
% mod_phase: (optional) matrix of size hslen*hrlen containing modulation 
%         phase values, which are supposed to tune the scale-rate filter 
%         to the signal (default is mod_factor=0)
% filt_type: string type, determines the moving direction of the filter 
%            'none' (full s-r domain)
%            'up' (upward analytic, nonzero over 2nd and 3rd quadrants)
%            'down' (downward analytic, nonzero over 1st and 4th quadrants)
%
% Outputs
% hsr: impulse response of the 2D filter (in time-frequency domain)
% Hsr: 2D-Fourier transform of the filter impulse response 
%
% Author: Fatemeh Pishdadian (fpishdadian@u.northwestern.edu)

%% extract filter parameters

% scale filter 
hs_len = scale_params.hslen; % length of the spectral filter
samprate_spec = scale_params.samprate_spec; % spectral ripple frequency
hs_type = scale_params.type; % type of the spectral filter

% rate filter 
beta = rate_params.time_const; % time constant of the temporal filter
hr_len = rate_params.hrlen;  % length of the temporal filter
samprate_temp = rate_params.samprate_temp; % # of frames per second
hr_type = rate_params.type;  % type of the temporal filter

%% frequency and time vectors

% zero-pad filters to the next even number 
hs_len = hs_len + mod(hs_len,2);
hr_len = hr_len + mod(hr_len,2);
    
% frequency and time vectors
freq_vec = (0:hs_len-1)'/samprate_spec;
time_vec = (0:hr_len-1)'/samprate_temp;

%% impulse response of the original scale filter: Gaussian 

hs = scale_ctr*(1-2*(scale_ctr*pi*freq_vec).^2).*exp(-((scale_ctr*freq_vec*pi).^2));
hs = [hs(1:hs_len/2+1);hs(hs_len/2:-1:2)]; % make it even so the transform is real

%% impulse response of the original rate filter    

hr = rate_ctr*(rate_ctr*time_vec).^2.*exp(-time_vec*beta*rate_ctr).*sin(2*pi*rate_ctr*time_vec);

%% scale response (Fourier transform of the scale impulse response)

% band-pass scale filter
Hs = abs(fft(hs,hs_len)); % discard the negligible imaginary part

% low/high-pass scale filter
if ~strcmp(hs_type,'bandpass')    
    Hs1 = Hs(1:hs_len/2+1,:); Hs1 = Hs1/max(Hs1);
    Hs2 = Hs(hs_len/2+2:end,:); Hs2 = Hs2/max(Hs2);
    [~,max_idx1] = max(Hs1);
    [~,max_idx2] = max(Hs2);
    
    if strcmp(hs_type,'lowpass')
        Hs1(1:max_idx1-1) = 1;  
        Hs2(max_idx2+1:end) = 1;
    elseif strcmp(hs_type,'highpass')
        Hs1(max_idx1+1:end) = 1;
        Hs2(1:max_idx2-1) = 1;
    end
    
    Hs = [Hs1;Hs2]; % form the full magnitude spectrum
end 

%% rate response (Fourier transform of the rate impulse response)

% band-pass rate filter
Hr = fft(hr,hr_len); % rate response is complex

% low/high-pass rate filter
if ~strcmp(hr_type,'bandpass')
    Hr_ph = angle(Hr); 
    Hr_mag = abs(Hr);
    
    Hr_mag1 = Hr_mag(1:hr_len/2+1,:); Hr_mag1 = Hr_mag1/max(Hr_mag1);
    Hr_mag2 = Hr_mag(hr_len/2+2:end,:); Hr_mag2 = Hr_mag2/max(Hr_mag2);
    [~,max_idx1] = max(Hr_mag1);
    [~,max_idx2] = max(Hr_mag2);
    
    if strcmp(hr_type,'lowpass')
       Hr_mag1(1:max_idx1-1) = 1;    
       Hr_mag2(max_idx2+1:end) = 1;
    elseif strcmp(hr_type,'highpass')
       Hr_mag1(max_idx1+1:end) = 1;
       Hr_mag2(1:max_idx2-1) = 1;
    end
    
    Hr_mag = [Hr_mag1;Hr_mag2]; % form the full magnitude spectrum
    Hr = Hr_mag.*exp(1j*Hr_ph); % form the full Fourier transfomr
end    
        
%% full scale-rate impulse and transform responses
 % Hsr_full is quadrant separable

Hsr_full = Hs*Hr.';

% normalize the filter magnitude:
Hsr_full_mag = abs(Hsr_full);
Hsr_full_mag = Hsr_full_mag/max(Hsr_full_mag(:));
Hsr_full = Hsr_full_mag.*exp(1j*angle(Hsr_full));

if strcmp(filt_dir,'up')
    
  % compute the upward version of the scale-rate response
  Hsr_up=Hsr_full;
  Hsr_up(1:hs_len/2,1:hr_len/2)=0; 
  %Hsr_up(hs_len/2+1:end,hr_len/2+1:end)=0; 
  Hsr_up(hs_len/2+2:end,hr_len/2+2:end)=0; 
  Hsr=Hsr_up;
  
elseif strcmp(filt_dir,'down') % downward ripple
    
  Hsr_down=Hsr_full;
%   Hsr_down(1:hs_len/2,hr_len/2+1:end)=0;
%   Hsr_down(hs_len/2+1:end,1:hr_len/2)=0;
  Hsr_down(1:hs_len/2,hr_len/2+2:end)=0;
  Hsr_down(hs_len/2+2:end,1:hr_len/2)=0;
  Hsr=Hsr_down;
  
else   
   Hsr=Hsr_full;
end

hsr=ifft2(Hsr);

if max(imag(hsr(:)))<1e-8, hsr=real(hsr); end

end


