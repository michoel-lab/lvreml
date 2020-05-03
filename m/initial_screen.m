function [beta2,varexpl,idx] = initial_screen(C,Z,varargin)
% INITIAL_SCREEN - Rapid screening of univariate variance component models
% INITIAL_SCREEN screens known covariates for possible filtering before
% applying the lvreml function. Screening is based on estimating the
% variance explained by each covariate alone. Details and rationale are in
% Section S4 of the paper.
% 
% USAGE: [beta2,varexpl] = initial_screen(C,Z)
%        [beta2,varexpl,idx] = initial_screen(C,Z,theta)
%
% INPUT: C - (n x n) empirical sample covariance matrix of expression data, 
%            see also data_prep
%        Z - (n x d) matrix of normalized data for d covariates (known
%            confounders) in n samples, see also data_prep
%
% OUTPUT: beta2   - (d x 1) vector, variance parameter in univariate lvreml
%                   model for each covariate, see Section S4 of the paper
%         varexpl - (d x 1) vector, variance in C explained by each
%                   covariate (= beta2/tr(C))
%         idx     - a linearly independent subset of covariates with 
%                   varexpl>theta, only returned if number of input arguments 
%                   equals 3  
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

[ns,nc] = size(Z); % number of covariates and samples
beta2 = zeros(nc,1);
varexpl = zeros(nc,1);
trC = sum(diag(C));
for k=1:nc
    beta2(k) = ns*Z(:,k)'*C*Z(:,k)/(ns-1) - trC/(ns-1);
end
% the solution only holds if beta2
varexpl(beta2>0) = beta2(beta2>0)./trC;


% get best linearly independent set with varexpl>threshold, if necessary
if nargout==3 && nargin==3
    vcut = varargin{1};
    [vs,t] = sort(varexpl,'descend');
    idx = t(vs>vcut);
    for k=2:length(idx)
        id = idx(1:k);
        rnk = rank(Z(:,id(id>0)));
        if rnk<sum(id>0)
            idx(k) = 0;
        end
    end
    idx = idx(idx>0);
end