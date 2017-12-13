function [X,alpha2,B,sigma2,nc2] = lvreml(Y,S,varexpl)
% LVREML - Restricted maximum-likelihood solution for linear mixed models with known and latent variance components 

[C,S] = data_prep(Y,S);
ns = length(C); % number of samples

% Full SVD of known covariates
[U,G,V] = svd(S);
nc2 = rank(G); % number of linearly independent covariates

% Split in space spanned by covariates and orthogonal complement
G1 = G(1:nc2,1:nc2);
U1 = U(:,1:nc2);
U2 = U(:,nc2+1:end);
V1 = V(:,1:nc2);

C11 = U1'*C*U1;
C22 = U2'*C*U2;

% Get reml estimate for X, pulled-back to original space
resvar = 1-varexpl;
if ~isempty(C22)
    [Vx,Ex] = eig(C22);
    [Evx,t] = sort(diag(Ex),'ascend');
    % variance not explained by X, at all possible numbers of hidden factors
    sigma2vec = cumsum(Evx)./(1:length(Evx))';
    nx = find(sigma2vec>resvar/(1-nc2*resvar/ns),1);
    X = U2*Vx(:,t(end:-1:nx));
    % Get estimate for residual and latent variable variances
    if nx>1
        sigma2 = mean(Evx(1:nx-1));
    else
        sigma2 = 0.;
    end
    alpha2 = Evx(end:-1:nx)-sigma2;
else
    X = [];
    alpha2 = [];
    sigma2 = 0;
end

% Get estimate for covariate covariance matrix
V1G = V1/G1; % V1*inv(G1);
B = V1G*(C11-sigma2*eye(nc2))*V1G';
