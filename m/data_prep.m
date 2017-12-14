function [C,Sn,Yn,idx] = data_prep(Y,S)
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

% Standardize Y and get overlap matrix
ng = size(Y,2); % number of genes
Yn = zscore(Y,1);
C = (Yn*Yn')/ng;

% Normalize covariates to have unit L2-norm
Sn = bsxfun(@rdivide,S,sum(S.^2).^.5);

% Check covariate matrix rank and return linearly independent subset if
% necessary
if nargout==4
    [~,R,idx] = qr(S,0);
    tol = max(size(S))*eps(max(abs(diag(R))));
    rk = sum(abs(diag(R))>tol); % rank
    if rk<size(S,2)
        idx = idx(1:rk);
    end
end