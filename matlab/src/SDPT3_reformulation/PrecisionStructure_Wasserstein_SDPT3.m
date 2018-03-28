function [ X, lambda, optimal_value ] = PrecisionStructure_Wasserstein_SDPT3( S_hat, rho, zero_mat, sqrt_S_hat )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Solve the precision matrix estimation using Wasserstein with given sparsity structure
%   The sparsity is given in the zero_mat: 
%   if zero_mat(i,j) > 0 then X(i,j) = 0
%

yalmip('clear');
p = size(S_hat, 1);

X_var = sdpvar(p, p);
A_var = sdpvar(p, p);

lambda_var = sdpvar(1);


CS = [lambda_var*eye(p) >= X_var >= 0, A_var >= 0, [A_var, lambda_var*sqrt_S_hat; lambda_var*sqrt_S_hat', lambda_var*eye(p) - X_var] >= 0];

zero_elements = find(zero_mat>0);
CS = [CS, X_var(zero_elements) == 0];

Obj = -logdet(X_var) + trace(A_var) + lambda_var*(rho^2 - trace(S_hat));


sdpt3options = sdpsettings('solver', 'sdpt3', 'debug', 0, 'verbose', 0);
diag_SDP = optimize(CS, Obj, sdpt3options);


% CS1 = [lambda_var*eye(p) - X_var >= 0, X_var >= 0];
% CS2 = [X_var >= 0];
% diag_SDP = optimize(CS2, trace(X_var*S_hat), sdpt3options)

X = value(X_var);
lambda = value(lambda_var);
optimal_value = value(Obj);
end

