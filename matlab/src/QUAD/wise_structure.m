function [ est ] = wise_structure( S_hat, rho, E, options )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
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
% 
% Output:
% est.value: the precision matrix estimate
% est.info: a termination message 

    if nargin == 4
        input.S = S_hat;
        input.rho = rho;
        input.E = E;
        
        est = wise_structure_main(input, options);
    else
        est = NaN;
        disp('Input not correct!');
    end
end

