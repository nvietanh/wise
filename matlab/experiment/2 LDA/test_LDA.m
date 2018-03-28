function [ nb_of_correct_prediction ] = test_LDA( mu, Precision, test_data )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% given the mu and Precision matrix, do the testing on the test data
% return the numbef of correct prediction
    
    % initialize
    nb_of_correct_prediction = 0;
    % get the number of class
    C = length(mu);
    
    for i = 1:length(test_data.y)
        score = zeros(C, 1);
        for c = 1:C
            score(c) = (test_data.X(i,:)-mu{c})*Precision*(test_data.X(i,:)-mu{c})';
        end
        [~, idx] = min(score);

        if abs(idx-test_data.y(i)-1)<0.1
            % correct prediction
            nb_of_correct_prediction = nb_of_correct_prediction + 1;
        end
    end
end

