function [ y ] = wise_cov_find_y(lambda, gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Covariance Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function solves a cubic equation to find y for the covariance
%   estimator


    w = nthroot((gamma*sqrt(lambda) + gamma*sqrt(lambda + 4*gamma/27))/2, 3);
    y = (w - gamma/(3*w))^2;

end

