function [ distance ] = WassersteinDistanceNormal( mu1, Sigma1, mu2, Sigma2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRO Precision Matrix Estimation
% Viet Anh NGUYEN, Peyman MOHAJERIN, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [V, D] = eig(Sigma1);
    sqrtSigma1 = V*sqrt(D)*V';
    [V, D] = eig(sqrtSigma1*Sigma2*sqrtSigma1);
    
    distance = norm(mu1-mu2, 2)^2 + trace(Sigma1) + trace(Sigma2) - 2*trace(sqrt(D));
    distance = sqrt(distance);

end

