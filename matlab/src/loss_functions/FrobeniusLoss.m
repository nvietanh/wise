function [ loss ] = FrobeniusLoss( Sigma1, Sigma2 )
%   The Frobenius norm between 2 covariance matrix, normalized by dimension
    loss = norm(Sigma1 - Sigma2, 'fro')/size(Sigma1,1);
    
end

