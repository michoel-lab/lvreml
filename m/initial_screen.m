function beta2 = initial_screen(Y,S)
% INITIAL_SCREEN - Rapid screening of univariate variance component models

[C,S] = data_prep(Y,S);
[ns,nc] = size(S); % number of covariates and samples
beta2 = zeros(nc,1);
for k=1:nc
    beta2(k) = (ns-2)*S(:,k)'*C*S(:,k)/(ns-1) - sum(diag(C))/(ns-1);
end