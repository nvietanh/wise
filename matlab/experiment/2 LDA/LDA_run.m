%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Linear Discriminant Analysis Application
% Leave m Out Cross Validation method

clear all; close all; clc;
% add the src path 
startup

%%
data_name = 'colon'; % 'colon' or 'leuk'

% Set m for the cross validation of the training data set
% if m = 1, we get the LOOCV
m = 1;
p_vector = [20, 40, 80, 100, 200];

% specify the number of samples used for training for each class
nb_training = [27, 11]; % for leuk data

if strcmp(data_name, 'colon')
    nb_training = [9, 20]; % for colon data
end

xfile = ['data/' data_name '_xtr.txt'];
yfile = ['data/' data_name '_ytr.txt'];
X_data = dlmread(xfile);
y_data = dlmread(yfile);
disp(['Dataset: ' data_name]);

C = length(unique(y_data)); % number of classes



%% Separate the data into a training set and a validation set
for i = 1:C
    class_idx{i} = find(y_data == i-1);
end


train_idx = [];
for i = 1:C
    train_idx = [train_idx; class_idx{i}(1:nb_training(i))];
end
test_idx = setdiff(1:length(y_data), train_idx);

train_X_data = X_data(train_idx, :);
train_y_data = y_data(train_idx);

valid_X_data = X_data(test_idx, :);
valid_y_data = y_data(test_idx);

%%
% run LmO-CV to learn the optimal parameters for each methods

rho_vect = 10.^[-1:0.05:2];
alpha_vect = 10.^[-3:0.05:0];
beta_vect = 10.^[-3:0.05:0];

n = length(train_y_data); %number of observation
correct_Wass = zeros(length(p_vector), length(rho_vect));
correct_Cov = zeros(length(p_vector), length(alpha_vect));
correct_l1 = zeros(length(p_vector), length(beta_vect));



% Divide the training data into n choose m parts

% Set of indices for the test set

test_idx = combnk(1:n, m);
train_idx = [];
% Set of indices for the training set
for k = 1:size(test_idx, 1)
    train_idx(k, :) = setdiff(1:n, test_idx(k,:));
end
    
for i = 1:length(p_vector)

    p = p_vector(i); % the number of genes used
      
    % Leave-2-Out Cross Validation
    for k = 1:size(test_idx, 1)
        % training set
        X = train_X_data(train_idx(k,:), 1:p);
        y = train_y_data(train_idx(k,:), :);
        % test set
        test_data.X = train_X_data(test_idx(k, :), 1:p);
        test_data.y = train_y_data(test_idx(k, :), :); 

        error = zeros(size(X));
        c_index{1} = find(y < 0.5); % index for class 0
        c_index{2} = find(y > 0.5); % index for class 1

        % calculate mean and error
        for c = 1:C
            mu{c} = mean(X(c_index{c}, :), 1);
            error(c_index{c}, :) = X(c_index{c}, :) - repmat(mu{c}, length(c_index{c}), 1);
        end
        % calculate the pooled sample covariance matrix
        S_hat = error'*error/(n-C); %the degree of freedom is n-C

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For Wasserstein robust
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:length(rho_vect)
            rho = rho_vect(j);
            % calculate the Wasserstein precision matrix estimator
            est = wise(S_hat, rho, 1e-8, 1e-5);
            Precision = est.value;
            % test performance on the test_data
            correct_Wass(i,j) = correct_Wass(i,j) + test_LDA(mu, Precision, test_data);
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For Shrinkage covariance estimation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:length(alpha_vect)
            alpha = alpha_vect(j);

            % Shrinkage estimation
            Precision = inv((1-alpha)*S_hat + alpha*diag(diag(S_hat)));
            correct_Cov(i,j) = correct_Cov(i,j) + test_LDA(mu, Precision, test_data);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For l_1 QUIC
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:length(beta_vect)
            beta = beta_vect(j);

            % Shrinkage estimation
            Precision = QUIC(S_hat, beta, 1e-6, 0);
            % correct prediction
            correct_l1(i,j) = correct_l1(i,j) + test_LDA(mu, Precision, test_data);
        end
       
    end % end of k Cross-Validation
    
    display(['Training iteration: ' num2str(i) ' over ' num2str(length(p_vector))]);
end % end of loop for the number of genes

% End of training sample
%%
% Do classification on the validation set
oos_correct_Wass = zeros(length(p_vector),1);
oos_correct_Shrink = zeros(length(p_vector),1);
oos_correct_l1 = zeros(length(p_vector),1);

for i = 1:length(p_vector)
    
    % First, calculate the best parameters
    [val, max_idx] = max(correct_Wass(i,:));
    rho = rho_vect(max_idx);

    [val, max_idx] = max(correct_Cov(i,:));
    alpha = alpha_vect(max_idx);

    [val, max_idx] = max(correct_l1(i,:));
    beta = beta_vect(max_idx);

    
    p = p_vector(i);
    % First, estimate mu and precision matrix from the test data
    % calculate mean and error
    X = train_X_data(:, 1:p);
    y = train_y_data;
    c_index{1} = find(y < 0.5); % index for class 1
    c_index{2} = find(y > 0.5); % index for class 2
    mu = [];
    error = [];
    for c = 1:C
        mu{c} = mean(X(c_index{c}, :), 1);
        error(c_index{c}, :) = X(c_index{c}, :) - repmat(mu{c}, length(c_index{c}), 1);
    end
    % calculate the pooled sample covariance matrix
    S_hat = error'*error/(n-C); %the degree of freedom is n-C

    % Find the precision matrix for each method
    est = wise(S_hat, rho, 1e-8, 1e-5);
    Precision_Wass = est.value;
    
    Precision_Shrink = inv((1-alpha)*S_hat + alpha*diag(diag(S_hat)));
    Precision_l1 = QUIC(S_hat, beta, 1e-6, 0);
    
    % create the validation data set
    validation.X = valid_X_data(:, 1:p);
    validation.y = valid_y_data;
    oos_correct_Wass(i) = test_LDA(mu, Precision_Wass, validation);
    oos_correct_Shrink(i) = test_LDA(mu, Precision_Shrink, validation);
    oos_correct_l1(i) = test_LDA(mu, Precision_l1, validation);
    
    display(['Out-of-sample iteration: ' num2str(i) ' over ' num2str(length(p_vector))]);
end

disp('End of Experiment')
%% Calculate final results

%save([data_name '_LDA.mat']);

oos_mis_Wass = length(validation.y) - oos_correct_Wass';
oos_mis_Shrink = length(validation.y) - oos_correct_Shrink';
oos_mis_l1 = length(validation.y) - oos_correct_l1';

% overall_mis = [oos_mis_Wass; oos_mis_Shrink; oos_mis_l1];
overall_correct_rate = [oos_correct_Wass'; oos_correct_Shrink'; oos_correct_l1']/length(validation.y)*100;