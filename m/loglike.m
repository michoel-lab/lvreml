function L = loglike(K,C)
% LOGLIKE - Log-likelihood for LVREML output

[V,D] = eig(K);
Ev = diag(D);

L = -sum(log(Ev)) - sum(diag(D.^(-1)*V'*C*V));
