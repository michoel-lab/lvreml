function L = loglike(K,C)
% LOGLIKE - Log-likelihood for LVREML output

[V,D] = eig(K);
Ev = diag(D);

L = -sum(log(Ev)) - sum(diag(inv(D)*V'*C*V));
%L = - sum(diag(inv(D)*V*C*V'));