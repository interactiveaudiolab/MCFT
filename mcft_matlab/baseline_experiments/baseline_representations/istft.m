function x_out = istft(X,winlen,wintype,ovp)

% This function computes the inverse STFT (inverse spectrogram) using
% the Overlap Addition (OLA) method
% 
% Inputs:
% X: one-sided STFT/spectrogram of the signal x 
% winlen: length of the window (default window type: Hamming)
% wintype: window type,string (rectangular,bartlett,hamming,hanning,blackman)
% ovp: overlap between adjacent windows in STFT analysis
%
% Ouputs:
% x_out: the reconstructed signal
   
%% Form the full spectrogram (-pi,pi]
    
  XX=conj(X(end-1:-1:2,:));
  X=[X;XX]; % 2-sided spectrogram
  Nt=size(X,2);
  hop=winlen-ovp;  % window hop size
  
  
%% Generate winow samples:

switch wintype
    case 'rectangular'
        W=ones(winlen,1);
    case 'bartlett'
        W=bartlett(winlen);
    case 'hamming'
        W=hamming(winlen,'periodic');
    case 'hann'
        W=hann(winlen,'periodic');
    case 'blackman'
        W=blackman(winlen);
end

 
%% Reconstruction through OLA:

%  x_out=zeros((Nt-1)*hop+winlen,1); 
%  norm_w=zeros((Nt-1)*hop+winlen,1);
%  or equivalently:    
 x_out=zeros(Nt*hop+ovp,1); 
 norm_w=zeros(Nt*hop+ovp,1);

 for i=1:Nt  
     range_temp=(i-1)*hop+1:(i-1)*hop+winlen;
     inv_x_temp=real(ifft(X(:,i))).*W;
     x_out(range_temp)=x_out(range_temp)+inv_x_temp(1:winlen);
     norm_w(range_temp)=norm_w(range_temp)+W.^2;
 end  

ind_large_vals = find(norm_w>eps);
x_out(ind_large_vals)=x_out(ind_large_vals)./norm_w(ind_large_vals);

% Take out the zero-padded sections in the beginning and
% end of the signal (zero-padded sections were added when
% computing the stft, see f_stft)

% if the reconstructed signal is still longer than the 
% original signal, it means it has zero padding at the end

if ovp>=hop
    
    ovp_hop_raio=ceil(ovp/hop);    
    z_start=ovp_hop_raio*hop;
    z_end=ovp;
    
    x_out(1:z_start)=[];
    x_out(end-z_end+1:end)=[];
    
elseif ovp<hop
    
    z_start=hop;    
    x_out(1:z_start)=[];
    
    z_end=ovp;
    x_out(end-z_end+1:end)=[];
    
            
end

 
end
  
         
         
         
         
         
         
         
         
