function [ Output ] = PlotTube(data, param, ColorTube)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot a tube for data given the size "Interval"
% data: Run * N

x_vect         = param.x_vect; % vector for values of rho
Interval  = param.Interval;

I = prctile(data,Interval,1);
Output.MIN  = I(1,:);
Output.MAX  = I(2,:);
Output.plot = fill([x_vect,fliplr(x_vect)], [Output.MIN fliplr(Output.MAX)], ColorTube);
set(Output.plot,'EdgeColor','None');

Output.AVR  = mean(data,1);

end

% data = matrix of dimensions (nb of runs,dim(rho))