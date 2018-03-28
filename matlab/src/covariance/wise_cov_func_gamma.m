function [ value ] = wise_cov_func_gamma( gamma, lambda, rho )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
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

