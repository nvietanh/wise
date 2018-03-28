function [ val ] = wise_structure_f( S_hat, rho, X, gamma, L )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate the true DRO objective function
% f(X, gamma) = -logdet(X) + gamma^2*Trace( (gamma*I - X)^{-1} S_hat ) + gamma*(rho^2 - Trace(S_hat))
%
% Input:
% S_hat: the sample covariance matrix
% rho: the Wasserstein ambiguity radius
% X: the precision matrix estimator
% gamma: the dual variable
% L: the Cholesky decomposition of X

    if nargin < 5
        L = chol(X);
    end
    
    p = size(S_hat, 1);
    
    val = - 2*sum(log(diag(L))) + gamma*(rho^2 - trace(S_hat) + trace((eye(p)-X/gamma)\S_hat));
end

