function [beta2,varexpl,idx] = initial_screen(Y,S,varargin)
% INITIAL_SCREEN - Rapid screening of univariate variance component models

[C,S] = data_prep(Y,S);
[ns,nc] = size(S); % number of covariates and samples
beta2 = zeros(nc,1);
trC = sum(diag(C));
for k=1:nc
    beta2(k) = ns*S(:,k)'*C*S(:,k)/(ns-1) - trC/(ns-1);
end
varexpl = beta2./trC;

% get best linearly independent set with varexpl>threshold, if necessary
if nargout==3 && nargin==3
    vcut = varargin{1};
    [vs,t] = sort(varexpl,'descend');
    idx = t(vs>vcut);
    for k=2:length(idx)
        id = idx(1:k);
        rnk = rank(S(:,id(id>0)));
        if rnk<sum(id>0)
            idx(k) = 0;
        end
    end
    idx = idx(idx>0);
end