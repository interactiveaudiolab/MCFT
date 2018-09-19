function X=inv_cft(cft_in,params)

% This function computes the inverse common fate transform of a 2d signal
%
% Inputs:
% cft_in: 4d matrix contatining the cft coefficients
% params: transform parameters including
%         1. spec_size (vector containing original spectrogram dimensions)
%         2. hop_size ((vector containing vertical and horizontal hop
%         sizes)
% 
% Output:
% X: 2d matrix containing the reconstructed spectrogram

%% Dimensions 
[Lw1,Lw2,Nv,Nh]=size(cft_in);

Nft=params.spec_size;
Hop=params.hop_size;

Win=[Lw1,Lw2];
Nvh=[Nv,Nh];
Ovp=Win-Hop;

Weights = ones(Lw1,Lw2);
%Weights = hamming(Lw1,'periodic')*hamming(Lw2,'periodic').';

Xz_dims=Nvh.*Hop+Ovp; % dimensions of the zero-padded spectrogram

%% Compute the inverse transform

Xz=zeros(Xz_dims(1),Xz_dims(2));
norm_w=zeros(size(Xz));
for i=1:Nv
    vrange=(i-1)*Hop(1)+(1:Win(1));
    for j=1:Nh
        hrange=(j-1)*Hop(2)+(1:Win(2));  
        Xz(vrange,hrange)=Xz(vrange,hrange)+ifft2(cft_in(:,:,i,j)).* Weights;
        norm_w(vrange,hrange)=norm_w(vrange,hrange)+ Weights.^2;
    end
end

X=Xz(Ovp(1)+1:Ovp(1)+Nft(1),Ovp(2)+1:Ovp(2)+Nft(2));
norm_w=norm_w(Ovp(1)+1:Ovp(1)+Nft(1),Ovp(2)+1:Ovp(2)+Nft(2));
X=X./norm_w;




