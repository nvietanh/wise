function [ alpha_best, newX, newGamma, new_obj, L ] = wise_structure_find_stepsize(input, sol, direction, delta, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Find the stepsize using Armijo's rule
%
%

    % initialize the parameters for line search
    beta = options.beta;
    alpha = options.alpha;
    sigma = options.sigma;
    max_linesearch_iter = options.max_linesearch_iter;
    
    found = 0;
    iter = 0;
    alpha_best = 0;
    new_obj = nan;
    while found < 0.5 && iter < max_linesearch_iter
        
        % Calculate the new solution
        newX = sol.X;
        newX(sol.idx) = newX(sol.idx) + alpha*direction.Delta_X;
        newX = triu(newX)+triu(newX,1)';
        
        newGamma = sol.gamma + alpha*direction.Delta_gamma;
        
        % do a cholesky on X
        [L, not_psd] = chol(newX, 'lower');
        if not_psd < 0.5 && eigs(newX, 1) < newGamma - options.eig_tol
            % if the new solution is feasible, then calculate the obj
            new_obj = wise_structure_f(input.S, input.rho, newX, newGamma, L);
            
            if new_obj < sol.f + alpha*sigma*delta 
                % if sufficient descent is met, then terminate
                found = 1;
                alpha_best = alpha;
            end
        end
        
        % decrease alpha
        alpha = beta*alpha;   
        
        iter = iter + 1;
    end

end

