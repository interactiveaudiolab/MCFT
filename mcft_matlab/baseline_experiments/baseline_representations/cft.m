function cft_out=cft(X,params)

% This function computes the common fate transform of a 2d signal
%
% Inputs:
% X: 2d matrix contatining the spectrogram of an audio signal 
% params: transform parameters including
%         1. win_size (vector containing height and width of the window)
%         2. hop_size (vector containing vertical and horizontal hop sizes)
% 
% Output:
% cft_out: 4d matrix containing cf coefficients

%% Dimensions and parameters
[Nf,Nt]=size(X);

Win=params.win_size;
Hop=params.hop_size;
Ovp=Win-Hop;

Weights=1;
%Weights = hamming(Win(1),'periodic')*hamming(Win(2),'periodic').';
%Weights = blackman(Win(1))*blackman(Win(2)).';


% compute the number of tiles 
Nvh=ceil(([Nf,Nt]-Ovp)./Hop);

% zero pad the signal
Nzft=(Nvh.*Hop)+3*Ovp; % size of the zero-padded signal
Xz=zeros(Nzft);
Xz(Ovp(1)+1:Ovp(1)+Nf,Ovp(2)+1:Ovp(2)+Nt)=X;

%% Compute the transform
Nzvh=floor((Nzft-Ovp)./Hop);

cft_out=zeros(Win(1),Win(2),Nzvh(1),Nzvh(2));
for i=1:Nzvh(1)
    vrange=(i-1)*Hop(1)+(1:Win(1));
    for j=1:Nzvh(2)
        hrange=(j-1)*Hop(2)+(1:Win(2));       
        X_temp=fft2(Xz(vrange,hrange).*Weights);
        cft_out(:,:,i,j) = X_temp;
    end
end


