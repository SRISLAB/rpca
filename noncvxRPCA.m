% This is the implementation of nonconvex RPCA method in ICDM 2015 paper
% 'Robust PCA via Nonconvex Rank Approximation'
% Zhao Kang, August 2015. Questions? Zhao.Kang@siu.edu;

clear all
clc

name='subject5';
load([name,'.mat']);
[m,n]=size(STFT);


muzero=.7;    % the only tuning parameter

lambda=4.5e-3;  % model parameter
type=21;   %different modeling of Sparsity.
rate=0.5;   %update rate of \mu
gamma=3e-3;     %gamma parameter in the rank approximation
tol=1e-2;  % stopping criterion

%initializations
S=zeros(m,n);
Y=zeros(m,n);
L=STFT;
sig=zeros(min(m,n),1); % for DC 
mu=muzero;

tic;
for ii=1:500
  
    D=STFT-S-Y/mu;
    [ L,sig] = DC(D,mu/2,sig,gamma);
    [S]=errorsol(Y,STFT,L,lambda,mu,type);
    Y=Y+mu*(L-STFT+S);
    mu=mu*rate;
    
    sigma=norm(STFT-S-L,'fro');
    RRE=sigma/norm(STFT,'fro');
    
    if RRE<tol
        break
    end
    
end
time_cost = toc;
rk=rank(L)

imagesc(abs(L))

