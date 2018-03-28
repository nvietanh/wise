function [ sol, msg ] = wise_structure_initialize( varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    if nargin == 2
        input = varargin{1};
        initial_sol = varargin{2};
    else
        msg = 'Error! Not enough input arguments to initialize solution!';
        disp(msg);
        return
    end
    
    % First, initialize input.E2
    % input.E2 keeps track of elements which can be non-zeros
    input.E2 = 1 - input.E;
    % E2 to take care of the symmetric elements (only take the upper part)
    input.E2(find(tril(input.E2,-1) > 0.1)) = 0; %reset lower part to 1


    % initialize solution
    sol.idx = find(input.E2>0.5);
    [sol.lat, sol.lon] = find(input.E2 > 0.5);
    sol.diag_element = +(sol.lat==sol.lon);
    sol.offdiag_element = +(sol.lat~=sol.lon);
    % p is the size of the covariance matrix
    sol.p = size(input.E2, 1); 
    % d is the number of nonzeros (diagonal + upper triangle) 
    % [effective size of the problem]
    sol.d = length(sol.idx); 

    sol.linear_weight = sol.diag_element + 2*sol.offdiag_element;

    if ~isstruct(initial_sol) || ~wise_structure_check_solution_consistency(input, initial_sol)
        % If there is no initial_sol, then initial X and gamma
        % Set X_0, gamma_0
        %sol.X = inv(diag(diag(input.S)));
        %sol.gamma = sqrt(trace(sol.X*input.S*sol.X))/input.rho;      
        sol.X = eye(sol.p);
        %sol.gamma = 1 + 1;
        sol.gamma = sqrt(trace(sol.X*input.S*sol.X))/input.rho;
    else
        % If initial_sol is provided, and that solution is consistent, initialize with the provided sol
        sol.X = initial_sol.X;
        sol.gamma = initial_sol.gamma;
       
    end
    
    % find X inverse
    sol.Xinv = inv(sol.X);
    
    sol.rho2_m_trs = input.rho^2 - trace(input.S);
    sol.f = wise_structure_f(input.S, input.rho, sol.X, sol.gamma, chol(sol.X, 'lower'));
        
    msg = 'Initialization successful!';
end

