function [ est ] = wise( varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the Wasserstein Inverse-covariance Shrinkage Estimator
%
    
    if nargin == 2
        S_hat = varargin{1};
        rho = varargin{2};
        eig_tol = 1e-14;            % default value
        bisection_tol = 1e-5;       % default value
    elseif nargin == 3
        S_hat = varargin{1};
        rho = varargin{2};
        eig_tol = varargin{3};
        bisection_tol = 1e-5;       % default value
    elseif nargin == 4
        S_hat = varargin{1};
        rho = varargin{2};
        eig_tol = varargin{3};
        bisection_tol = varargin{4};
    else
        disp('Error! Input arguments mismatch');
        return
    end
    
    [W, D] = mexeig_dgesdd(S_hat);
   
    lambda = diag(D);
    p = length(lambda);
    
    good_idx = find(lambda > eig_tol);
    null_idx = find(lambda <= eig_tol);
    
    if length(good_idx) < p && rho < eig_tol
        % Singular sample covariance matrix and rho close to 0
        est.value = NaN;
        est.msg = 'Singular sample covariance matrix, rho too close to 0. No precision matrix estimator available';
        est.gamma = NaN;
        est.eig = NaN;
        est.eigbasis = NaN;
        est.min_gamma = NaN;
        est.max_gamma = NaN;
    else
        lambda_bar = lambda(good_idx);
        lambda_max = max(lambda);

        % Find the range for lambda search
        min_gamma = max(p*(p*lambda_max + 2*rho^2 - p*lambda_max*sqrt(1 + 4*rho^2/p/lambda_max))/(2*rho^4), 1e-14);
             
        if length(good_idx) < p
            max_gamma = p/rho^2;
        else
            max_gamma = min(p/rho^2, sqrt(sum(1./lambda_bar))/rho);
        end
        
        % When rho is too small, it may happen that min_lambda > max_lambda
        % Correct if this happens
        if min_gamma > max_gamma
            min_gamma = eig_tol;
        end
        
        % Calculate optimal gamma from bisection
        gamma = wise_bisect(@wise_func_gamma, lambda_bar, rho, [min_gamma, max_gamma], bisection_tol);

        % Calculate the eigenvalues
        x = zeros(p, 1);     % eigenvalues of the precision matrix estimator
        for i = 1:length(good_idx)
            h = gamma*lambda(good_idx(i));
            
            % The below expression seems to be the most stable
            new_val = gamma*(1 - 2/(1+sqrt(1+4/h)));
            
            % If the above calculation is not stable, then try a new
            % expression
            if new_val <= 0
                % recalculate new_val
                new_val = gamma*(1 + (h*(1 - sqrt(1 + 4/h)))/2);
            end
            
            if new_val <= 0
                % recalculate new_val
                sqrt_h = sqrt(h);
                new_val = gamma*(1 + sqrt_h*(sqrt_h - sqrt(h+4))/2);
            end
            
            if new_val <= 0
                new_val = eig_tol;
                disp('Unstable computation of eigenvalues!');
            end
            
            
            x(good_idx(i)) = new_val;
        end

        for i = 1:length(null_idx)
            x(null_idx(i)) = gamma;
        end
        
        est.value = W*diag(x)*W';
        est.eig = x;
        est.eigbasis = W;
        est.gamma = gamma;
        est.min_gamma = min_gamma;
        est.max_gamma = max_gamma;
        est.msg = 'Precision matrix estimator available';
    end
end
