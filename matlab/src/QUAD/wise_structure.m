function [ est ] = wise_structure( varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the Wasserstein Inverse-covariance Shrinkage Estimator with
% sparsity constraints
%
% Input:
% S_hat: sample covariance matrix
% rho: Wasserstein ambiguity set radius
% E: matrix specifying the sparsity, E(i,j) > 0 means X(i,j) = 0. The upper
% triangular of E is used. 
% options: an options structure
% X_init: initial gamma for X (optional)
% gamma_init: initial value for gamma (optional)
% 
% Output:
% est.value: the precision matrix estimate
% est.info: a termination message 

    if nargin == 4
        input.S = varargin{1};
        input.rho = varargin{2};
        input.E = varargin{3};
        options = varargin{4};
        est = wise_structure_main(input, options);
    elseif nargin == 6
        input.S = varargin{1};
        input.rho = varargin{2};
        input.E = varargin{3};
        options = varargin{4};
        initial_sol.X = varargin{5};
        initial_sol.gamma = varargin{6};
        est = wise_structure_main(input, options, initial_sol);
    else
        est = NaN;
        disp('Input not correct!');
    end
end

