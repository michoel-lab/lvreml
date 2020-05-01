function [C,Zn,Yn] = data_prep(Y,Z)
% DATA_PREP - Prepare expression and known covariate data for use in lvreml functions
% DATA_PREP prepares the expression and known covariate data for use in
% other functions of the lvreml package. Details and rationale are in
% Section S2 of the paper.
% 
% USAGE: [C,Zn,Yn] = data_prep(Y,Z)
%
% INPUT: Y - (n x m) matrix of expression data for m genes in n samples
%        Z - (n x d) matrix of data for d covariates (known confounders) in
%            n samples
%
% OUTPUT: C  - empirical sample covariance matrix for the expression data,
%         Zn - normalized covariate data (all variables scaled to have unit
%              L2 norm)
%         Yn - centred expression data (all samples centred to have mean zero)   
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

% Check that data is in right format (columns are variables, rows are samples)
if ~isempty(Z)
    if size(Y,1)~=size(Z,1)
        % rows are not matching, check columms
        if size(Y,2)==size(Z,2)
            % colums are samples -> transpose both
            Y = Y';
            Z = Z';
        else
            % inconsistent data format
            error('lvreml::data_prep::cannot determine sample dimension');
        end
    end
else
   warning('lvreml::data_prep::no covariates provided, assuming expression data is in format samples x genes'); 
end

% Center samples (rows of Y) to remove fixed effects on mean
Yn = bsxfun(@minus,Y,mean(Y,2));
% Compute overlap matrix
ng = size(Y,2); % number of genes
C = (Yn*Yn')/ng;

% Normalize covariates to have unit L2-norm
if ~isempty(Z)
    Zn = bsxfun(@rdivide,Z,sum(Z.^2).^.5);
else
    Zn = [];
end