function [L,S] = adm(M,u,lambda,itemax,tol)
% solve PCP problem by ADM algorithm
% v1.0 2020-1-1
% function：solve the following optimization problem
%                  min  ||X||*+lambda||E||_F
%                  s.t. M = A+E

% initialize S0 and Y0 and L0
[m,n] = size(M) ;
S = zeros(m,n) ;
Y = S ;
L = S ;

% the observed matrix can contain non number
unobserved = isnan(M);
M(unobserved) = 0;

if nargin < 2
    lambda = 1 / sqrt(max(m,n));
end
if nargin < 3
    u = 10*lambda;
end
if nargin < 4
    tol = 1e-6;
end
if nargin < 5
    itemax = 1000;
end

for ii = 1:itemax
    L = sig_thre(M-S+(1/u)*Y,(1/u)) ;
    S = soft_thre(M-L+(1/u)*Y, lambda/u) ;
    Z = M-L-S ;
    Y = Y+u*Z ;
    error = norm(M-L-S,'fro')/norm(M,'fro') ;
    if (ii == 1) || (mod(ii, 10) == 0) || (error < tol)
        fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
            ii, error, rank(L), nnz(S));
    end
    if error<tol
        break;
    end
end

