function [ y ] = wise_cov_find_y(lambda, gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This function solves a cubic equation to find y for the covariance
%   estimator


    w = nthroot((gamma*sqrt(lambda) + gamma*sqrt(lambda + 4*gamma/27))/2, 3);
    y = (w - gamma/(3*w))^2;

end

