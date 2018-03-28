%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This script run the synthetic simulation for the simulation result in the
% paper
%
% In the paper we use the density number d, but in this script, we use the
% variables nb_nonzeros to change the density
%
% nb_nonzeros = 50  corresponds to d = 12.5%
% nb_nonzeros = 100 corresponds to d = 25%
% nb_nonzeros = 200 corresponds to d = 50%
% nb_nonzeros = 400 corresponds to d = 100%

%%
clear all; close all; clc;
startup;

%%
% Generate matrix
rng(0,'v5uniform');
p = 20;
nb_nonzeros = 100; % refer to the description to set this value



U = zeros(p);
nonzeros_idx = randperm(p^2, nb_nonzeros);
U(nonzeros_idx) = sign(randn(nb_nonzeros,1));

A = U'*U;

eigA = eig(A);
Sigma_inv = A + max(1e-4, -1.2*min(eigA))*eye(p);
Sigma = inv(Sigma_inv);

zi = find(abs(triu(Sigma_inv))<1e-6);
[lat lon] = find(abs(Sigma_inv)<1e-6);
good_idx = find(lat< lon);
lat = lat(good_idx);
lon = lon(good_idx);

% specify the matrix of zero elements
pct_info = [0.5, 0.75, 1];
for k = 1:length(pct_info)
    to_keep = ceil(pct_info(k)*length(lat));
    lat1 = lat(1:to_keep);
    lon1 = lon(1:to_keep);
    zero_mat = zeros(p);
    for i = 1:length(lat1)
        zero_mat(lat1(i), lon1(i)) = 1;
        zero_mat(lon1(i), lat1(i)) = 1;
    end
    param.zero_mat{k} = zero_mat;
end

NbInstance = 100;
seed_vect = randi(1e6, [max(1000,NbInstance), 1]);
 
param.DataSampleSize_vect = [10, 20, 40, 60];

% for the wasserstein estimator
param.rho_vect = 10.^[-5:0.05:1];
% for the QUIC
param.beta_vect = 10.^[-4:0.04:0];
% for the linear shrinkage estimator
param.alpha_vect = 1e-5:1e-5:1;

param.Sigma = Sigma;

% Define the loss function
param.LossFunction = @SteinLoss;

%% full run
for instance_count = 1:NbInstance
    param.SeedNumber = seed_vect(instance_count);

    [Error] = synthetic_simulation(param);
    AllError.Wasserstein_Robust_Error(:, :, instance_count) = Error.Wasserstein_Robust_Error;
    AllError.l1_Error_SDPT3(:, :, instance_count) = Error.l1_Error_SDPT3;
    AllError.Covariance_Shrinkage_Error(:, :, instance_count) = Error.Covariance_Shrinkage_Error;
    AllError.Wasserstein_Structure_Error(:, :, :, instance_count) = Error.Wasserstein_Structure_Error;
    display(['Current progress: ' num2str(instance_count) ' over ' num2str(NbInstance)]);
end



