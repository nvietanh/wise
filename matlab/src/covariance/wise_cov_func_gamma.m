function [ value ] = wise_cov_func_gamma( gamma, lambda, rho )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Covariance Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function used for bisection (used to find lambda in the covariance matrix problem)
%
    value = rho^2 - sum(lambda(:));
    for i = 1:length(lambda)
        y = wise_cov_find_y(lambda(i), gamma);
        value = value + y*(1+2*y/gamma);
    end

end

