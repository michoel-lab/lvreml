function [X,alpha2,B,D,sigma2,K] = lvreml(Y,Z,targetX)
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

    % Get REML estimate for X
    if ~isempty(C22)
        [Vx,Ex] = eig(C22);
        [Evx,t] = sort(diag(Ex),'ascend');
        % minimum eigenvalue on known covariate space, needed to check validity
        % of the solution
        lambdamin = min(eig(C11));
        % variance not explained by X, at all possible numbers of latent
        % variables; this is called "f(p)" in the appendix of the paper
        sigma2vec = cumsum(Evx)./(1:length(Evx))';
        
        % check if 3rd arguments wants a given number of latent variables or
        % target variance explained; in both cases adjust to make sure the
        % solution will be valid
        if targetX >= 1
            nx = round(targetX); % number of latent variables
            % find the smallest number of latent variables we must include
            % (largest number we can cut)
            ncut_max = find(sigma2vec<lambdamin & Evx>Evx(1),1,'last');
            % cut fewer if we've asked for more latent variables
            ncut = min(length(Evx)-nx,ncut_max);
        elseif targetX>= 0
            varexpl = targetX; % target variance explained
            % Set the target residual variance
            resvar = min((1-varexpl)*trC/ns,lambdamin);
            % variance not explained by X, at all possible numbers of latent
            % variables; this is called "f(p)" in the appendix of the paper
            ncut = find(sigma2vec<resvar & Evx>Evx(1),1,'last');
        else
            error('lvreml::lvreml::3rd argument must be integer or value between zero and one');
        end
        if isempty(ncut)
            ncut = 0;
            sigma2 = 0.;
        else
            sigma2 = mean(Evx(1:ncut));
        end
        % pull back selected eigenvectors of C22 to original space
        X = U2*Vx(:,t(end:-1:ncut+1));
        alpha2 = Evx(end:-1:ncut+1)-sigma2;
    else
        X = [];
        alpha2 = [];
        sigma2 = 0;
    end
    
    % Get estimates for the covariance matrix among the known confounders
    V1G = V1/G1; % V1*inv(G1);
    B = V1G*(C11-sigma2*eye(nc))*V1G';
    
    % Get estimates for the covariance matrix between the known confounders and
    % latent variables
    if ~isempty(X)
        D = V1G*U1'*C*X;
    else
        D = [];
    end
    % Return the covariance matrix K
    K = [Z X]*[B D;D' diag(alpha2)]*[Z'; X'] + sigma2*eye(ns);
    K = 0.5*(K+K');
else
   % no covariates provided, only estimating hidden components using PCA
   B = [];
   D = [];
   [Vx,Ex] = eig(C);
   [Evx,t] = sort(diag(Ex),'ascend');
   Vx = Vx(:,t);
   % variance not explained by X, at all possible numbers of latent
   % variables; this is called "f(p)" in the appendix of the paper
   sigma2vec = cumsum(Evx)./(1:length(Evx))';
   % check how many we need to include at least
   nmin = find(Evx(2:end)-sigma2vec(1:end-1)>0,1,'first');
   if targetX >= 1
       nx = round(targetX); % number of latent variables
   elseif targetX>= 0
       varexpl = targetX; % target variance explained
       % set the target residual variance
       resvar = (1-varexpl)*trC/ns;
       % find required number of latent variables
       nx = find(sigma2vec<resvar,1,'last');
   else
       error('lvreml::lvreml::3rd argument must be integer or value between zero and one');
   end
   nx = max(nx,nmin);
   sigma2 = mean(Evx(1:end-nx));
   X = Vx(:,end:-1:end-nx+1);
   alpha2 = Evx(end:-1:end-nx+1)-sigma2;
   % Return the covariance matrix K
   K = X*diag(alpha2)*X' + sigma2*eye(ns);
   K = 0.5*(K+K');
end
