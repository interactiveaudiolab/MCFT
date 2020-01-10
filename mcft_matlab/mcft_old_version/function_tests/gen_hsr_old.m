function [hsr,Hsr]=gen_hsr_old(s,r,s_params,r_params,filt_type)

% This function generates a 2d rate/scale impulse response 
% 
% Inputs:
% s: scale value
% r: rate value 
% s_params:structure array containing thr paameters of the scale filte
%          hslen: length of the filter impulse response
%          ripple_freq: ripple frequency
%          type: string argument indicating the filter type
%               ('bandpass','lowpass','highpass')
% r_params:structure array containing thr paameters of the rate filte
%          time_const: exponent coefficient of the exponential term
%          hrlen: length of the rate filter impulse response
%          frame_per_sec: # of frames per second
%          type: string argument indicating the type of the filter
%                ('bandpass','lowpass','highpass')
%          r_support: string type, indicates the support of the 
%               filter in the rate domain. the filter can be nonzero 
%               over only positive rates, or only negative reates.
%               ('neg','pos')
% mod_phase: (optional) matrix of size hslen*hrlen containing modulation 
%         phase values, which are supposed to tune the scale-rate filter 
%         to the signal (default is mod_factor=0)
% filt_type: (optional) string type, determines the type of filter as 
%            'full' (over s-r domain) or 'analytic' (default is analytic)
%
% 
% Outputs
% hsr: impulse response of the scale-rate filter
% Hsr: Fourier transform of the filter impulse response (full spectrum)

%% set the filter type to 'analytic' if it is not provided
if nargin<5
    filt_type='analytic';
end

%% extract filter parameters

% scale filter 
Ls=s_params.hslen;
SRF=s_params.ripple_freq;
s_type=s_params.type;

% rate filter 
beta=r_params.time_const;
Lr=r_params.hrlen;
FPS=r_params.frame_per_sec;
r_type=r_params.type;
r_rng=r_params.r_support;

%% frequency and time vectors

% zero-pad filters to the nextpow2 
% Ls=2^nextpow2(Ls);
% Lr=2^nextpow2(Lr);
Ls=Ls+mod(Ls,2);
Lr=Lr+mod(Lr,2);
    
% frequency and time vectors
w = (0:Ls-1)'/SRF;
t = (0:Lr-1)'/FPS;

%% impulse response of the original scale filter: Gaussian 

hs=s*(1-2*(pi*s*w).^2).*exp(-(pi*s*w).^2);
%s=sqrt(2)*pi*s;
%hs=s*(1-(s*w).^2).*exp(-((s*w).^2)/2);

%hs=s*(s*w).^2.*exp(-w*beta*s).*sin(2*pi*s*w);
%hs=hs-mean(hs);

hs=[hs(1:Ls/2+1);hs(Ls/2:-1:2)]; % make it even so the transform is real

%% impulse response of the original rate filter    

hr=r*(r*t).^2.*exp(-t*beta*r).*sin(2*pi*r*t);
%hr=hr-mean(hr);

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

Hsr_full=Hs*Hr.';

% normalize the filter magnitude:
Hsr_full_mag=abs(Hsr_full);
Hsr_full_mag=Hsr_full_mag/max(Hsr_full_mag(:));
Hsr_full=Hsr_full_mag.*exp(1j*angle(Hsr_full));

if strcmp(filt_type,'analytic')
    
   % compute the analytic version of the scale-rate response
   Hsr_analytic=Hsr_full;
   %Hsr_analytic(Ls+1:end,:)=0;  %%% putting this to zero=imperfect reconst

   if strcmp(r_rng,'neg') % upward ripple
      Hsr_analytic(1:Ls/2,1:Lr/2)=0; 
      Hsr_analytic(Ls/2+2:end,Lr/2+2:end)=0; 
   elseif strcmp(r_rng,'pos') % downward ripple
      Hsr_analytic(1:Ls/2,Lr/2+2:end)=0;
      Hsr_analytic(Ls/2+2:end,1:Lr/2)=0;
   end
   
   hsr_analytic=ifft2(Hsr_analytic);
   hsr=hsr_analytic;
   if max(imag(hsr(:)))<1e-8, hsr=real(hsr); end
   Hsr=Hsr_analytic;
   
else
   hsr=hsr_full;
   if max(imag(hsr(:)))<1e-8, hsr=real(hsr); end
   Hsr=Hsr_full;
end

end


