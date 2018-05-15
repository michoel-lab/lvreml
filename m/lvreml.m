function [X,alpha2,B,sigma2] = lvreml(Y,Z,varexpl)
% LVREML - Restricted maximum-likelihood solution for linear mixed models with known and latent variance components 

[C,Z] = data_prep(Y,Z);
ns = length(C); % number of samples
trC = sum(diag(C));

% Full SVD of known covariates
if ~isempty(Z)
    [U,G,V] = svd(Z);
    nc = rank(G); % number of linearly independent covariates
    if nc~=size(Z,2)
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

% Set the target residual variance
resvar = min((1-varexpl)*trC/ns,min(eig(C11)));

% Get REML estimate for X
if ~isempty(C22)
    [Vx,Ex] = eig(C22);
    [Evx,t] = sort(diag(Ex),'ascend');
    % variance not explained by X, at all possible numbers of latent
    % variables; this is called "f(p)" in the appendix of the paper
    sigma2vec = cumsum(Evx)./(1:length(Evx))';
    nx = find(sigma2vec<resvar & Evx>Evx(1),1,'last');
    if isempty(nx)
        nx = 0;
        sigma2 = 0.;
    else
        sigma2 = mean(Evx(1:nx));
    end
    % pull back selected eigenvectors of C22 to original space
    X = U2*Vx(:,t(end:-1:nx+1));
    alpha2 = Evx(end:-1:nx+1)-sigma2;
else
    X = [];
    alpha2 = [];
    sigma2 = 0;
end

% Get estimate for the covariance matrix of the known confounders
if ~isempty(Z)
    V1G = V1/G1; % V1*inv(G1);
    B = V1G*(C11-sigma2*eye(nc))*V1G';
else
    B = [];
end
