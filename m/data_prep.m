function [C,Zn,Yn] = data_prep(Y,Z)
% DATA_PREP - Get expression data sample overlap and normalize covariates

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

% Center variables (columns of Y) to remove fixed effects on mean
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