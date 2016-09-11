function [hsr,Hsr]=gen_hsr(S,R,S_params,R_params,filt_type)

% This function generates a 2D-impulse response in the time-frequncy
% domain with dilation factors S and R: h(omega,tau;S,R)
%
% Inputs:
% S: scale value
% R: rate value 
% S_params:structure array containing the parameters of the spectral filter
%          hslen: length of the spectral filter impulse response
%          ripple_freq: ripple frequency (cyc/oct)
%          type: string argument indicating the filter type
%               ('bandpass','lowpass','highpass')
% R_params:structure array containing the parameters of the temporal filter
%          time_const: exponent coefficient of the exponential term
%          hrlen: length of the temporal filter impulse response
%          frame_per_sec: # of frames per second
%          type: string argument indicating the type of the filter
%                ('bandpass','lowpass','highpass')
% mod_phase: (optional) matrix of size hslen*hrlen containing modulation 
%         phase values, which are supposed to tune the scale-rate filter 
%         to the signal (default is mod_factor=0)
% filt_type: string type, determines the type of filter as 
%            'full' (full s-r domain)
%            'up' (upward analytic, nonzero over 2nd and 3rd quadrants)
%            'down' (downward analytic, nonzero over 1st and 4th quadrants)
%
% Outputs
% hsr: impulse response of the 2D filter (in time-frequency domain)
% Hsr: 2D-Fourier transform of the filter impulse response 

%% extract filter parameters

% scale filter 
Ls=S_params.hslen; % length of the spectral filter
SRF=S_params.ripple_freq; % spectral ripple frequency
s_type=S_params.type; % type of the spectral filter

% rate filter 
beta=R_params.time_const; % time constant of the temporal filter
Lr=R_params.hrlen;  % length of the temporal filter
FPS=R_params.frame_per_sec; % # of frames per second
r_type=R_params.type;  % type of the temporal filter

%% frequency and time vectors

% zero-pad filters to the next even number 
Ls=Ls+mod(Ls,2);
Lr=Lr+mod(Lr,2);
    
% frequency and time vectors
w = (0:Ls-1)'/SRF;
t = (0:Lr-1)'/FPS;

%% impulse response of the original scale filter: Gaussian 

hs=S*(1-2*(S*pi*w).^2).*exp(-((S*w*pi).^2));
hs=[hs(1:Ls/2+1);hs(Ls/2:-1:2)]; % make it even so the transform is real

%% impulse response of the original rate filter    

hr=R*(R*t).^2.*exp(-t*beta*R).*sin(2*pi*R*t);

%% scale response (Fourier transform of the scale impulse response)

% band-pass scale filter
Hs=abs(fft(hs,Ls)); % discard the negligible imaginary part

% low/high-pass scale filter
if ~strcmp(s_type,'bandpass')    
    Hs1=Hs(1:Ls/2+1,:); Hs1=Hs1/max(Hs1);
    Hs2=Hs(Ls/2+2:end,:); Hs2=Hs2/max(Hs2);
    [~,max_idx1]=max(Hs1);
    [~,max_idx2]=max(Hs2);
    
    if strcmp(s_type,'lowpass')
        Hs1(1:max_idx1-1)=1;  
        Hs2(max_idx2+1:end)=1;
    elseif strcmp(s_type,'highpass')
        Hs1(max_idx1+1:end)=1;
        Hs2(1:max_idx2-1)=1;
    end
    
    Hs=[Hs1;Hs2]; % form the full magnitude spectrum
end 

%% rate response (Fourier transform of the rate impulse response)

% band-pass rate filter
Hr=fft(hr,Lr); % rate response is complex

% low/high-pass rate filter
if ~strcmp(r_type,'bandpass')
    Hr_ph=angle(Hr); 
    Hr_mag=abs(Hr);
    
    Hr_mag1=Hr_mag(1:Lr/2+1,:); Hr_mag1=Hr_mag1/max(Hr_mag1);
    Hr_mag2=Hr_mag(Lr/2+2:end,:); Hr_mag2=Hr_mag2/max(Hr_mag2);
    [~,max_idx1]=max(Hr_mag1);
    [~,max_idx2]=max(Hr_mag2);
    
    if strcmp(r_type,'lowpass')
       Hr_mag1(1:max_idx1-1)=1;    
       Hr_mag2(max_idx2+1:end)=1;
    elseif strcmp(r_type,'highpass')
       Hr_mag1(max_idx1+1:end)=1;
       Hr_mag2(1:max_idx2-1)=1;
    end
    
    Hr_mag=[Hr_mag1;Hr_mag2]; % form the full magnitude spectrum
    Hr=Hr_mag.*exp(1j*Hr_ph); % form the full Fourier transfomr
end    
        
%% full scale-rate impulse and transform responses
 % Hsr_full is quadrant separable

Hsr_full=Hs*Hr.';

% normalize the filter magnitude:
Hsr_full_mag=abs(Hsr_full);
Hsr_full_mag=Hsr_full_mag/max(Hsr_full_mag(:));
Hsr_full=Hsr_full_mag.*exp(1j*angle(Hsr_full));

if strcmp(filt_type,'up')
    
  % compute the upward version of the scale-rate response
  Hsr_up=Hsr_full;
  Hsr_up(1:Ls/2,1:Lr/2)=0; 
  Hsr_up(Ls/2+2:end,Lr/2+2:end)=0; 
  Hsr=Hsr_up;
  
elseif strcmp(filt_type,'down') % downward ripple
    
  Hsr_down=Hsr_full;
  Hsr_down(1:Ls/2,Lr/2+2:end)=0;
  Hsr_down(Ls/2+2:end,1:Lr/2)=0;
  Hsr=Hsr_down;
  
else   
   Hsr=Hsr_full;
end

hsr=ifft2(Hsr);
if max(imag(hsr(:)))<1e-8, hsr=real(hsr); end

end


