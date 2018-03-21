function [ X, optimal_value ] = PrecisionSelection_SDPT3( S_hat, beta )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Solve the ell_1 regularized precision matrix estimation using SDPT3
%

    yalmip('clear');
    p = size(S_hat, 1);

    X_var = sdpvar(p, p);

    CS = [X_var >= 0];
    Obj = -logdet(X_var) + trace(S_hat*X_var) + beta*sum(sum(abs(X_var)));

    sdpt3options = sdpsettings('solver', 'sdpt3', 'debug', 0, 'verbose', 0);
    diag_SDP = optimize(CS, Obj, sdpt3options);

    X = value(X_var);
    optimal_value = value(Obj);
end

