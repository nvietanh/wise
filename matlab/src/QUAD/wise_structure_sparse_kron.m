function [ val ] = wise_structure_sparse_kron( A, B, offdiag_element, lon, lat )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function implement A kronecker B
% made specific for this problem to take into account the sparsity of X
%


    off = repmat(offdiag_element, 1, length(offdiag_element));
    
    val = A(lon, lon).*B(lat, lat);
    val = val + (A(lon, lat).*off').*(B(lat, lon).*off');
    val = val + (off.*A(lat, lon)).*(off.*B(lon, lat));
    val = val + (off.*A(lat, lat).*off').*(off.*B(lon, lon).*off');
end

