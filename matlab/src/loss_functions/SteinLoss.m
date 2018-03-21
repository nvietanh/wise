function [ loss ] = SteinLoss( Sigma, X )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Return the Stein's Loss as defined in the paper
% loss = -log(det(X*Sigma)) + trace(Sigma*X) - size(Sigma, 1);
    loss = - 2*sum(log(diag(chol(X, 'lower')))) - 2*sum(log(diag(chol(Sigma, 'lower')))) + trace(Sigma*X) - size(Sigma, 1);


end

