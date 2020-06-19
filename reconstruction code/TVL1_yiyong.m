function [ Recon1_tvspr ] = TVL1_yiyong(A_mat,PPa1,TVWeight,xfmWeight,Itnlim,noneg)
%%reference: Y. Han et.al, "Sparsity-based acoustic inversion in cross-sectional multiscale optoacoustic imaging," Medical physics 42, 5444-5452 (2015)

n=sqrt(size(A_mat,2));
xsize=[n,n];
W = @(x) Wavedb1(x,xsize,0); 
WT = @(x) Wavedb1(x,xsize,1);
TV = @(x)D(x);
GTV = @(x)adjD(x);

% initialize Parameters for reconstruction
param = init;
param.N=n;
param.FT = A_mat;
param.XFM = WT;
param.XFMIN = W;
param.TV = TV;
param.GTV =GTV;
param.data = PPa1;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
param.noneg=noneg; % 1 for nonnegative reconstruction, 0 for normal reconstruction


res=zeros(n,n);
res=W(res);

%res = TVL1_v1(res,param);
res = TVL1_LS_BB(res,param);
Recon1_tvspr = WT(res);

end

