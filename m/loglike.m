function L = loglike(K,C)
% LOGLIKE - Log-likelihood for sample covariance matrix
% LOGLIKE computes the (scaled) log-likelihood value for a sample covariance
% estimate K and the empirical sample covariance matrix C. See eq. (2) of
% the paper. 
% 
% USAGE: L = loglike(K,C)
%
% INPUT: K - (n x n) estimate for a sample covariance matrix
%        C - (n x n) empirical sample covariance matrix
%
% OUTPUT: L - log-likelihood value
%
% AUTHOR: Tom Michoel
%         tom.michoel@uib.no
%         https://lab.michoel.info
%
% REFERENCE: MA Malik and T Michoel. Restricted maximum-likelihood method
% for learning latent variance components in gene expression data with
% known and unknown confounders. 
%
% LICENSE: GNU GPL v3

[V,D] = eig(K);
Ev = diag(D);

L = -sum(log(Ev)) - sum(diag(inv(D)*V'*C*V));
