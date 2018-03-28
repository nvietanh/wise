function [ direction, delta, H, g] = wise_structure_find_direction( input, sol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% True DRO objective function
% Find the descent direction by solving a quadratic approx problem
% Already with reduced direction to reduce memory
% This version is implementing x'*H*x + g*x (not 0.5*x'*H*x + g*x as in
% paper)

    
    % Contrary to the paper, in this code, G is defined as the inverse already
    G = inv(eye(sol.p) - sol.X/sol.gamma);
    GSG = G*input.S*G;
    GSGXG = GSG*sol.X*G;
    H_aux = -(G*sol.X*GSG + GSGXG)/sol.gamma^2;
    
    % Construct H and g
    H = zeros(sol.d+1);
    H(1:sol.d, 1:sol.d) = wise_structure_sparse_kron(sol.Xinv, sol.Xinv, sol.offdiag_element, sol.lon, sol.lat)/2 + 1/sol.gamma*wise_structure_sparse_kron(GSG, G, sol.offdiag_element, sol.lon, sol.lat);
    H(1:sol.d, sol.d+1) = H_aux(sol.idx).*sol.linear_weight/2;
    H(sol.d+1, 1:sol.d) = H(1:sol.d, sol.d+1)';
    H(sol.d+1, sol.d+1) = trace(GSGXG*sol.X)/sol.gamma^3;
    
      
    g(1:sol.d,1) = (GSG(sol.idx)-sol.Xinv(sol.idx)).*(sol.linear_weight);
    g(sol.d+1,1) = (sol.rho2_m_trs +trace(G*input.S) - trace(GSG*sol.X)/sol.gamma);

    % One can apply a block matrix inversion here, but it does not help much
    % A = H(1:sol.d, 1:sol.d);
    % B = H(1:sol.d, sol.d+1);
    % C = B';
    % D = H(end,end);
    % invH = [ inv(A-B*C/D), -inv(A-B*C/D)*B/D; -C*inv(A-B*C/D)/D, 1/D+C*inv(A-B*C/D)*B/D^2];
    % optimal_y = -invH*g/2;
    
    optimal_y = -H\g/2;
    direction.Delta_X = optimal_y(1:(end-1));
    direction.Delta_gamma = optimal_y(end);
    
    % compute delta for this direction
    delta = optimal_y'*g;
end

