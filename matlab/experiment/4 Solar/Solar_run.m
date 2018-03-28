%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The solar irradiation application
%
clear all; close all; clc;
startup;

%%
load(['data' filesep 'diablerets_data_kfold.mat']);

%%
options = wise_structure_settings('iter_limit', 100, 'sigma', 1e-4, 'delta_tol', 1e-6, 'gradient_tol', 1e-4);
%
rho_vect_Wstruct =  10.^(linspace(-2, 0, 60));
Stein_Wstruct = nan(length(rho_vect_Wstruct), length(kfold));
Fro_Wstruct = nan(length(rho_vect_Wstruct), length(kfold));
FroPrec_Wstruct = nan(length(rho_vect_Wstruct), length(kfold));
for i = 1:length(rho_vect_Wstruct)
    display(['rho ' num2str(i)]);
    rho = rho_vect_Wstruct(i);
    for k = 1:length(kfold)
        S_hat = kfold{k}.training;
        est = wise_structure(S_hat, rho, E, options);
                
        Stein_Wstruct(i,k) = SteinLoss(kfold{k}.test, est.value);
        Fro_Wstruct(i,k) = FrobeniusLoss(kfold{k}.test, inv(est.value));         % Frobenius in the covariance space
        FroPrec_Wstruct(i,k) = FrobeniusLoss(est.value, inv(kfold{k}.test));     % Frobenius in the precision space
    end
end

%%
% Covariance Shrinkage estimation
alpha_vect = 0.0001:0.0001:0.03;
Stein_SAA = nan(length(alpha_vect), length(kfold));
Fro_SAA = nan(length(alpha_vect), length(kfold));
FroPrec_SAA = nan(length(alpha_vect), length(kfold));
for i = 1:length(alpha_vect)
    display(['alpha ' num2str(i)]);
    alpha = alpha_vect(i);
    for k = 1:length(kfold)
            S_aux = (1-alpha)*kfold{k}.training + alpha*diag(diag(kfold{k}.training));
            X = inv(S_aux);
            Stein_SAA(i,k) =  SteinLoss(kfold{k}.test, X);
            Fro_SAA(i,k) = FrobeniusLoss(kfold{k}.test, S_aux);
            FroPrec_SAA(i,k) = FrobeniusLoss(X, inv(kfold{k}.test));
    end
end
%%
% l1 estimation
beta_vect = 10.^[-5:0.1:-3];
Stein_QUIC = nan(length(beta_vect), length(kfold));
Fro_QUIC = nan(length(beta_vect), length(kfold));
FroPrec_QUIC = nan(length(beta_vect), length(kfold));
% after the first run with beta_vect as above, we restrict to a smaller
% region as below:
%beta_vect = 10.^[-4.1:0.01:-3.9];
for i = 1:length(beta_vect)
    display(['beta ' num2str(i)]);
    beta = beta_vect(i);
    for k = 1:length(kfold)
        X = QUIC(kfold{k}.training, beta, 1e-6, 0);
        Stein_QUIC(i,k) = SteinLoss(kfold{k}.test, X);
        Fro_QUIC(i,k) = FrobeniusLoss(kfold{k}.test, inv(X));
        FroPrec_QUIC(i,k) = FrobeniusLoss(X, inv(kfold{k}.test));
    end
end

%%
% Wasserstein estimator without structure
rho_vect_W =  10.^[-2:0.01:0];
Stein_W = nan(length(rho_vect_W), length(kfold));
Fro_W = nan(length(rho_vect_W), length(kfold));
FroPrec_W = nan(length(rho_vect_W), length(kfold));

for i = 1:length(rho_vect_W)
    display(['rho ' num2str(i)]);
    rho = rho_vect_W(i);
    for k = 1:length(kfold)
        S_hat = kfold{k}.training;
        
        est = wise(S_hat, rho, 1e-7, 1e-5);
        X = est.value;
        
        Stein_W(i,k) = SteinLoss(kfold{k}.test, X);
        Fro_W(i,k) = FrobeniusLoss(kfold{k}.test, inv(X));         % Frobenius in the covariance space
        FroPrec_W(i,k) = FrobeniusLoss(X, inv(kfold{k}.test));     % Frobenius in the precision space
    end    
end
%%
display('done');

