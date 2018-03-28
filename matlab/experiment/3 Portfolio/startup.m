%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add the necessary path
path1 = pwd;
path2 = strsplit(path1,filesep);
path3 = strjoin(path2(1:end-2),filesep);
path_src = strcat(path3,filesep,'src');
path(genpath(path_src),path); 
clear path1 path2 path3 path_src