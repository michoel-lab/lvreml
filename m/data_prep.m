function [C,Sn] = data_prep(Y,S)
% DATA_PREP - Get expression data sample overlap and normalize covariates

% Check that data is in right format (columns are variables, rows are samples)
if size(Y,1)~=size(S,1)
    % rows are not matching, check columms
    if size(Y,2)==size(S,2)
        % colums are samples -> transpose both
        Y = Y';
        S = S';
    else
        % inconsistent data format
        error('lvreml::data_prep::cannot determine sample dimension');
    end
end

% Center Y and get overlap matrix
ng = size(Y,2); % number of genes
Y = bsxfun(@minus,Y,mean(Y));
C = (Y*Y')/ng;

% Normalize covariates to have unit L2-norm
Sn = bsxfun(@rdivide,S,sum(S.^2).^.5);