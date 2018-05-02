function [X,alpha2,B,sigma2,per,loglike] = lvreml(Y,S,varexpl)
% LVREML - Restricted maximum-likelihood solution for linear mixed models with known and latent variance components 

[C,S,~,per] = data_prep(Y,S);
ns = length(C); % number of samples

% Full SVD of known covariates
if ~isempty(S)
    [U,G,V] = svd(S(:,per));
    nc = rank(G); % number of linearly independent covariates
    if nc~=length(per)
        error('lvreml::lvreml::inconsitent rank estimates');
    end
    
    % Split in space spanned by covariates and orthogonal complement
    G1 = G(1:nc,1:nc);
    U1 = U(:,1:nc);
    U2 = U(:,nc+1:end);
    V1 = V(:,1:nc);
    
    C11 = U1'*C*U1;
    C22 = U2'*C*U2;
else
   % no covariates provided, only estimating hidden components
   nc = 0;
   C11 = [];
   C22 = C;
   U2 = eye(size(C));
end

% Log-likelihood for space spanned by covariates
E1 = eig(C11);
loglike = [-sum(log(E1))-nc, 0];

% Get reml estimate for X, pulled-back to original space
resvar = 1-varexpl;
if ~isempty(C22)
    [Vx,Ex] = eig(C22);
    [Evx,t] = sort(diag(Ex),'ascend');
    % variance not explained by X, at all possible numbers of hidden factors
    sigma2vec = cumsum(Evx)./sum(Evx);%./(1:length(Evx))';
    nx = find(sigma2vec>resvar,1);
    X = U2*Vx(:,t(end:-1:nx));
    % Get estimate for residual and latent variable variances
    if nx>1
        sigma2 = mean(Evx(1:nx-1));
        lsig = log(sigma2);
    else
        sigma2 = 0.;
        lsig = 0.;
    end
    alpha2 = Evx(end:-1:nx)-sigma2;
    % Log-likelihood for residual space
    loglike(2) = -sum(log(Evx(t(end:-1:nx)))) - (nx-1)*lsig - ns + nc;
else
    X = [];
    alpha2 = [];
    sigma2 = 0;
end

% Get estimate for covariate covariance matrix
if ~isempty(S)
    V1G = V1/G1; % V1*inv(G1);
    B = V1G*(C11-sigma2*eye(nc))*V1G';
else
    B = [];
end
