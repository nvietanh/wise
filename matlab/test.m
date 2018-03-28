%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate data
rng(2018, 'v5uniform');

% Set dimension p and number of observations n
p = 20;
n = 15;

% Generate the true covariance matrix
Sigma = randn(p);
Sigma = Sigma*Sigma'; L = chol(Sigma);
% Generate sample data using standard Gaussian and calculate sample covariance
data = L*randn(p, n);
S = data*data'/n;

% Ambiguity radius
rho = 0.5;
%% 
% Wasserstein Inverse-covariance Shrinkage Estimate
est1 = wise(S, rho);

eigen_tolerance = 1e-14;
bisection_tolerance = 1e-1;

est2 = wise(S, rho, eigen_tolerance, bisection_tolerance);

% Find the error, notice that the true covariance matrix is the identity
error1 = SteinLoss(Sigma, est1.value)
error2 = SteinLoss(Sigma, est2.value)


%%

% Set up the sparsity information
E = zeros(p);

% Only need to specify the sparsity in the upper triangle of E
E(1, 2) = 1;

options = wise_structure_settings('iter_limit', 500, 'sigma', 1e-4, 'delta_tol', 1e-6, 'gradient_tol', 1e-5, 'verbose', 1);
est_struct = wise_structure(S, rho, E, options);
error_struct = SteinLoss(Sigma, est_struct.value)