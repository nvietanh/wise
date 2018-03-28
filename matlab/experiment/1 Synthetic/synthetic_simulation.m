function [ Error ] = synthetic_simulation(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Sigma = param.Sigma;
    SeedNumber = param.SeedNumber;
    DataSampleSize_vect = param.DataSampleSize_vect;
    rho_vect = param.rho_vect; 
    beta_vect = param.beta_vect;
    alpha_vect = param.alpha_vect;


    % find the dimension for the problem
    p = size(Sigma,1);

    % set seed and generate data
    rng(SeedNumber,'v5uniform')
    data_all = mvnrnd(zeros(p,1), Sigma, max(max(DataSampleSize_vect), 1000));

    Error.Wasserstein_Robust_Error = zeros(length(DataSampleSize_vect), length(rho_vect));
    %Error.SA_Error = zeros(length(DataSampleSize_vect), 1);
    Error.Covariance_Shrinkage_Error = zeros(length(DataSampleSize_vect), length(alpha_vect));
    Error.Wasserstein_Structure_Error = zeros(length(param.zero_mat), length(DataSampleSize_vect), length(rho_vect));
    Error.l1_Error_SDPT3 = zeros(length(DataSampleSize_vect), length(beta_vect));

    for i = 1:length(DataSampleSize_vect)
        data = data_all(1:DataSampleSize_vect(i),:);
        % Compute the empirical covariance
        S_hat = data'*data/size(data,1);
        
        % Calculate error for sample average is not stable 
        %Error.SA_Error(i) = param.LossFunction(Sigma, inv(S_hat));

        % Run Wasserstein DRO estimation
        for idx_rho = 1:length(rho_vect)
            rho = rho_vect(idx_rho);
            est = wise(S_hat, rho, 1e-14, 1e-5);
            Error.Wasserstein_Robust_Error(i, idx_rho) = param.LossFunction(Sigma, est.value);   
        end

        % Run l1-estimation - introduced by Aspremont
        for idx_beta = 1:length(beta_vect)
            beta = beta_vect(idx_beta);
            X = QUIC(S_hat, beta, 1e-6, 0);
            Error.l1_Error(i, idx_beta) = param.LossFunction(Sigma, X);
        end

        % Run covariance shrinkage-estimator    
        Target = diag(diag(S_hat));
        for idx_alpha = 1:length(alpha_vect)
            alpha = alpha_vect(idx_alpha);
            Est = (1-alpha)*S_hat + alpha*Target;
            Error.Covariance_Shrinkage_Error(i, idx_alpha) = param.LossFunction(Sigma, inv(Est));
        end


        % run Wasserstein with structure
        [V, D] = eig(S_hat);
        D = sqrt(max(D, 0));
        sqrt_S_hat = V*D*V';
        sqrt_S_hat = (sqrt_S_hat+sqrt_S_hat')/2; %make it symmetric
        % Run SDPT3 to find the estimator
        for k = 1:length(param.zero_mat)
            for idx_rho = 1:length(rho_vect)
                rho = rho_vect(idx_rho);
                [ X  ] = PrecisionStructure_Wasserstein_SDPT3(S_hat, rho, param.zero_mat{k}, sqrt_S_hat);
                Error.Wasserstein_Structure_Error(k, i, idx_rho) = param.LossFunction(Sigma, X);
            end
        end


        %solve l_1 using SDPT3 (not necessarily the same output as QUIC)
        %Run SDPT3 to find the estimator
        for idx_beta = 1:length(beta_vect)
              beta = beta_vect(idx_beta);
              [ X  ] = PrecisionSelection_SDPT3(S_hat, beta);
              Error.l1_Error_SDPT3(i, idx_beta) = param.LossFunction(Sigma, X);
        end

    end % end of for loop


end

