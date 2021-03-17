function [X,alpha2,B,D,sigma2,K] = lvreml(C,Z,targetX)
% LVREML - Restricted maximum-likelihood solution for linear mixed models with known and latent variance components 
% LVREML computes the restricted maximum-likelihood solution for linear mixed 
% models with known and latent variance components  Details and rationale are in
% Sections S2, S4-S7 of the paper.
%
% IMPORTANT NOTE: Input parameters C and Z must be scaled and normalized
% using the "data_prep" function to ensure correct output!
% 
% USAGE: [X,alpha2,B,D,sigma2,K] = lvreml(C,Z,targetX)
%
% INPUT: C       - (n x n) sample covariance matrix of expression data for
%                  m genes in n samples, see also data_prep
%        Z       - (n x d) matrix of L1-normalized data for d covariates
%                  (known confounders) in n samples, see also data_prep
%        targetX - number, interpreted as target number of latent
%                  variables (targetX >= 1), or target value for the total
%                  variance in Y explained by the model (0<=targetX<1)
%
% OUTPUT: X      - (n x p) matrix of latent variable data (rows are
%                  samples, columns are variables, p is automatically 
%                  determined), also called X in the paper.
%         alpha2 - (p x 1) vector of variance explained by each latent
%                  variable; this is the diagonal of the diagonal matrix A
%                  in the paper.  
%         B      - (d x d) matrix of covariances among the known
%                  covariates, also called B in the paper.
%         D      - (d x p) matrix of covariances between known and latent
%                  variables, also called D in the paper.
%         sigma2 - the squared residual variance
%         K      - the restricted maximum-likelihood sample covariance
%                  matrix
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

ns = length(C); % number of samples
trC = sum(diag(C));

% Full SVD of known covariates
if ~isempty(Z)
    [U,G,V] = svd(Z);
    nc = rank(G); % number of linearly independent covariates
    if nc~=size(Z,2)
        error('lvreml::lvreml::inconsistent rank estimates');
    end
    
    % Split in space spanned by covariates and orthogonal complement
    G1 = G(1:nc,1:nc);
    U1 = U(:,1:nc);
    U2 = U(:,nc+1:end);
    V1 = V(:,1:nc);
    
    % Corresponding matrix blocks
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
        
        % check if 3rd argument wants a given number of latent variables or
        % target variance explained; in both cases adjust to make sure the
        % solution will be valid
        if targetX >= 1
            nx = round(targetX); % number of latent variables
            % find the smallest number of latent variables we must include
            % (largest number we can cut)
            ncut_max = find(sigma2vec<lambdamin & Evx>Evx(1),1,'last');
            % cut fewer if we've asked for more latent variables (but never
            % less than 0)
            ncut = min(max(length(Evx)-nx,0),ncut_max);
        elseif targetX>= 0
            varexpl = targetX; % target variance explained
            % Set the target residual variance
            resvar = min((1-varexpl)*trC/ns,lambdamin);
            % find number of variables that need to be cut to explain no
            % less than the target variance explained
            ncut = find(sigma2vec<resvar & Evx>Evx(1),1,'last');
        else
            error('lvreml::lvreml::3rd argument must be integer or value between zero and one');
        end
        if isempty(ncut)
            ncut = 0;
            sigma2 = 0.;
        else
            % estimate sigma^2 parameter
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
    
    % Get estimates for the covariance matrix between the known confounder
    % and latent variable effects
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
   % check how many we need to include at least (largest possible value of
   % excluded dimensions)
   nmax = find(Evx(2:end)-sigma2vec(1:end-1)>0,1,'last');
   % check how many we want to include from user input
   % note that nx is the number of dimensions *not* to include
   if targetX >= 1
       nx = ns - round(targetX); % "nx = n-p"
   elseif targetX>= 0
       varexpl = targetX; % target variance explained
       % set the target residual variance
       resvar = (1-varexpl)*trC/ns;
       % find required number of latent variables
       nx = find(sigma2vec<resvar,1,'last'); % "nx = n-p"
   else
       error('lvreml::lvreml::3rd argument must be integer or value between zero and one');
   end
   nx = min(nx,nmax);
   sigma2 = mean(Evx(1:nx));
   X = Vx(:,end:-1:nx+1);
   alpha2 = Evx(end:-1:nx+1)-sigma2;
   % Return the covariance matrix K
   K = X*diag(alpha2)*X' + sigma2*eye(ns);
   K = 0.5*(K+K');
end
