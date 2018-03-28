function [ est ] = wise_structure_main( varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The main function performing the sequential quadratic programming
%   problem for inverse covariance estimation.
%
%   input.E contains the location of zeros in the precision matrix
%
    if nargin == 2
        input = varargin{1};
        options = varargin{2};
        initial_sol = NaN;
    elseif nargin == 3
        input = varargin{1};
        options = varargin{2};
        initial_sol = varargin{3};
    else
        disp('Error! Input arguments mismatch');
        return
    end
    
    tic
    display(['Optimization using quadratic approximation.']);
    if options.verbose
        display(['----------------------------------------------']);
        display(sprintf('%6s\t%12s\t%12s\t%12s\t%10s','Iter','Stepsize', 'NormGrad','Obj.Value','Time'));
    end
    
    [sol, msg] = wise_structure_initialize(input, initial_sol);
    
    if options.verbose
        display(sprintf('%6i\t%10e\t%15s\t%10e\t%10.2f', 0, 0, '-    ', sol.f, toc));
    end
    
    iter = 1;
    
    while iter > 0
        
        % Find direction
        [direction, delta] = wise_structure_find_direction(input, sol);
        % Find step size
        [alpha, newX, newGamma, new_obj] = wise_structure_find_stepsize(input, sol, direction, delta, options);
        % Update solution
        if(~isnan(new_obj))
            sol = wise_structure_update_solution(sol, newX, newGamma, new_obj);
        end
            
        info.obj(iter) = sol.f;
        info.stepsize(iter) = alpha;
        info.norm_grad(iter) = norm(direction.Delta_X, 2)/sol.d + norm(direction.Delta_gamma);
        info.toc(iter) = toc;
        % Display iteration information
        if options.verbose
            wise_structure_print_output(info);
        end
        
        
        % Check termination criteria
        if info.norm_grad(iter) < options.gradient_tol 
            sol.message = 'Gradient norm close to 0. Optimal solution found.';
            break;
        elseif abs(delta)/sol.d < options.delta_tol
            sol.message = 'delta close to 0. Optimal solution found.';
            break;
        elseif toc > options.time_limit 
            sol.message = 'Time limit reached.'; 
            break;
        elseif iter > options.iter_limit-1
            sol.message = 'Iteration limit reached.';
            break;
        elseif isnan(new_obj)
            sol.message = 'Step size too small.';
            break;
        end    

        
        iter = iter+1;
    end
    display(sol.message);
    disp('Optimization procedure terminated.');     
    
    est.value = sol.X;
    est.info  = sol.message;
end

