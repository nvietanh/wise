function [ consistent ] = wise_structure_check_solution_consistency( input, initial_sol )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   This function is used to check whether the initial solution provided by
%   the user is consistent with the feasible set

    consistent = true;
    
    % Check symmetry and zero elements
    if issymmetric(initial_sol.X) || norm(input.E(:)'*initial_sol.X(:), 2) > 1e-10
        consistent = false;
    end
    
    % check lambda*I succeq X succeq 0
    [L, not_psd] = chol(initial_sol.X, 'lower');
    if not_psd > 0.5 && eigs(initial_sol.X, 1) < initial_sol.lambda
        consistent = false;
    end
    
end

