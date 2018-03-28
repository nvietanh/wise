function [ sol ] = wise_structure_update_solution(sol, newX, newGamma, new_obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Update solution to newX, newGamma
%

    % update
    sol.X = newX;
    sol.gamma = newGamma;
    sol.f = new_obj;
    sol.Xinv = inv(sol.X);   
end

