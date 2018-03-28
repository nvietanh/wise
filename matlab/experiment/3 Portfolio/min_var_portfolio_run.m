%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Portfolio Optimization Application
%

clear all; close all; clc;
startup


%%
dataset = '100'; % '100' or '48'
% Number of observations for the training data set
% also equal to the number of historical data to estimate the precision
% matrix
n = 120;          


%load data
load(['data' filesep dataset '.mat']);
scale_factor = 100;
return_data = return_data/scale_factor;

end_training_date = 200512;           
% Period for balancing
balancing_period = 3;
% Get the number of assets
p = size(return_data,2);

%%
% Run LOOCV on the training data set

end_train_idx = find(time_data == end_training_date);

rho_vect = 10.^[-2:0.01:0];
alpha_vect = 10.^[-2:0.01:0];
beta_vect = 10.^[-4:0.02:0];
o = ones(p, 1);

W_return = zeros(n, length(rho_vect));
S_return = zeros(n, length(alpha_vect));
Q_return = zeros(n, length(beta_vect));

idx = (end_train_idx-n+1):end_train_idx;

train_data = return_data(idx, :);
idx = 1:n;
for k = 1:n
    
    train_idx = find(idx ~= k);
    test_idx = k;
    
    % calculate S_hat
    data = train_data(train_idx, :);
    mu = mean(data,1);
    centered_data = data - repmat(mu, size(data,1), 1);
    S_hat = centered_data'*centered_data/(n-1);
    
    % Wasserstein precision
    for i = 1:length(rho_vect)
        rho = rho_vect(i);
        est = wise(S_hat, rho, 1e-10, 1e-4);
        X = est.value;
        w = X*o/sum(sum(X));
        W_return(k, i) = train_data(k,:)*w;
        disp('here')
    end
    
    % Covariance shrinkage
    for i = 1:length(alpha_vect)
        alpha = alpha_vect(i);
        invX = (1-alpha)*S_hat + alpha*diag(diag(S_hat));
        w = invX\o/(o'*(invX\o));
        S_return(k,i) = train_data(k,:)*w;
        disp('there')
    end
    
    
    % ell_1 regularization
    for i = 1:length(beta_vect)
        beta = beta_vect(i);
        X = QUIC(S_hat, beta, 1e-6, 0);
        w = X*o/sum(sum(X));
        Q_return(k,i) = train_data(k,:)*w;
        disp('bere')
    end
    
    
    display(['Fold: ' num2str(k) ' of ' num2str(n)]);
end

%% Find optimal value for alpha, beta, rho

mean_W = mean(W_return, 1);
mean_S = mean(S_return, 1);
mean_Q = mean(Q_return, 1);

centered_W = W_return - repmat(mean_W, size(W_return, 1), 1);
centered_S = S_return - repmat(mean_S, size(S_return, 1), 1);
centered_Q = Q_return - repmat(mean_Q, size(Q_return, 1), 1);

std_W = sqrt(diag(centered_W'*centered_W/(size(W_return, 1)-1)))';
std_S = sqrt(diag(centered_S'*centered_S/(size(S_return, 1)-1)))';
std_Q = sqrt(diag(centered_Q'*centered_Q/(size(Q_return, 1)-1)))';

[val, idx_W] = min(std_W);
[val, idx_S] = min(std_S);
[val, idx_Q] = min(std_Q);

rho = rho_vect(idx_W);
alpha = alpha_vect(idx_S);
beta = beta_vect(idx_Q);


%%
% The financial crisis ends on 200906, which is time 162
start_time = end_train_idx+1;
end_time = size(return_data, 1);



Wasserstein_return = zeros(size(return_data,1),1);
Shrinkage_return = zeros(size(return_data,1), 1);
QUIC_return = zeros(size(return_data,1), 1);

for t = start_time:end_time
    if mod(t, balancing_period) == 1 || balancing_period < 2
        %need to re-estimate the mean and precision matrix
        display('Rebalancing');
        data = return_data((t-n):(t-1), :);
        mu = mean(data,1);
        centered_data = data - repmat(mu, size(data,1), 1);
        S_hat = centered_data'*centered_data/(n-1);
        
        % re-estimate X
        est = wise(S_hat, rho, 1e-10, 1e-6);
        X = est.value;
        w_wasserstein =X*o/sum(sum(X));
        
        % re-estimate X
        invX = (1-alpha)*S_hat + alpha*diag(diag(S_hat));
        w_shrink = invX\o/(o'*(invX\o));
        
        % re-estimate X
        X = QUIC(S_hat, beta, 1e-6, 0);
        w_quic = X*o/sum(sum(X));
    end

    Wasserstein_return(t) = return_data(t,:)*w_wasserstein;
    Shrinkage_return(t) = return_data(t,:)*w_shrink;   
    QUIC_return(t) = return_data(t,:)*w_quic;
    display(['Out-of-sample iteration: ' num2str(t)]);
end
% End of computation


%% Result analysis
final_W = Wasserstein_return(start_time:end_time)*scale_factor;
final_S = Shrinkage_return(start_time:end_time)*scale_factor;
final_Q = QUIC_return(start_time:end_time)*scale_factor;

disp('Mean return: ');
mean(final_W)
mean(final_S)
mean(final_Q)

disp('Standard deviation: ');
std(final_W)
std(final_S)
std(final_Q)

disp('Sharpe ratio: ');
mean(final_W)/std(final_W)
mean(final_S)/std(final_S)
mean(final_Q)/std(final_Q)


% stack into one matrix
res(1,:) = [std(final_W), mean(final_W), mean(final_W)/std(final_W), rho];
res(2,:) = [std(final_S), mean(final_S), mean(final_S)/std(final_S), alpha];
res(3,:) = [std(final_Q), mean(final_Q), mean(final_Q)/std(final_Q), beta];

%%
% save(['FF' dataset '_n' num2str(n) '.mat']);