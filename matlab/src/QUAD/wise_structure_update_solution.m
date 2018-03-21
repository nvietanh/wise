function [ sol ] = wise_structure_update_solution(sol, newX, newGamma, new_obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
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

