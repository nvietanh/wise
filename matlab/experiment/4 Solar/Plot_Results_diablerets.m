%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is to plot the graph used in the paper
% Run this script after running the SolarPrecision.m script
% One can also download the result file from the first author's website


%%
% Load the result file (if available)
% Download the result file from the first author's website and put the file
% into the same folder


% clear all; close all; clc;
% load('result/diablerets_result.mat');


%% Calculate the average values
average_Stein_QUIC    = mean(Stein_QUIC, 2);
average_Stein_SAA     = mean(Stein_SAA, 2);
average_Stein_W       = mean(Stein_W, 2);
average_Stein_Wstruct = mean(Stein_Wstruct, 2);
%%
% plotting parameters
FontSize  = 24; FontSize2  = 18;
LineWidth = 2;  LineWidth2 = .1;
MarkerSize = 8;

Trans = [0.15, 0.25, 0.35, 0.45, 0.55];

param_plot.Interval = [0,  100];

cc=[0.0000    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4660    0.6740    0.1880
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
  
ColorTube(1,:)  = cc(1,:);
ColorTube(2,:)  = cc(2,:);
ColorTube(3,:)  = cc(3,:);
ColorTube(4,:)  = cc(4,:);
ColorTube(5,:)  = cc(5,:);

LineMarkType = {'-.', ':k', '--', ':', '.'};

%%

rho_toplot_Wstruct   = 1:length(rho_vect_Wstruct)-5;    % Wasserstein with structure
rho_toplot_W  = 120:length(rho_vect_W)-20;              % Wasserstein without structure

figure;
hold on; 
combined_legends = [];

count = 1;
% first, plot Wasserstein without structure
param_plot.x_vect       = rho_vect_W(rho_toplot_W);
plot(param_plot.x_vect, average_Stein_W(rho_toplot_W), LineMarkType{count}, 'Color', ColorTube(count,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
combined_legends{count} = ['Wasserstein shrinkage w/o structure'];

count = 4;
% forth, plot Wasserstein with structure
param_plot.x_vect       = rho_vect_Wstruct(rho_toplot_Wstruct);
plot(param_plot.x_vect, average_Stein_Wstruct(rho_toplot_Wstruct), LineMarkType{count}, 'Color', ColorTube(count,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
combined_legends{2} = ['Wasserstein shrinkage w structure'];

% now, plot the tube
% plot for Wasserstein without structure
toplot = 1;
param_plot.x_vect       = rho_vect_W(rho_toplot_W);
P_GP   = PlotTube(Stein_W(rho_toplot_W, :)', param_plot, ColorTube(toplot,:));
set(P_GP.plot,'facealpha',Trans(toplot));

% plot for Wasserstein with structure, true objective
toplot = 4;
param_plot.x_vect       = rho_vect_Wstruct(rho_toplot_Wstruct);
P_GP   = PlotTube(Stein_Wstruct(rho_toplot_Wstruct, :)', param_plot, ColorTube(toplot,:));
set(P_GP.plot,'facealpha',Trans(toplot));

% add legend
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', 15, 'Location', 'NorthWest')

ylim([250, 500]);
ylabel('Stein''s loss','Interpreter','latex','FontSize',FontSize2);
xlabel('$\rho$','Interpreter','latex','FontSize',FontSize2);
%title('Performance of 3 types of estimator','Interpreter','latex','FontSize',FontSize2);
set(gca,'ygrid','on')
%set(gca,'yscale','log')
set(gca,'xscale','log')

% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%%
alpha_toplot = 1:length(alpha_vect); % shrinkage

figure;
hold on; 
combined_legends = [];
count = 2;
% second, plot Covariance Shrinkage
param_plot.x_vect       = alpha_vect(alpha_toplot);
plot(param_plot.x_vect, average_Stein_SAA(alpha_toplot), LineMarkType{count}, 'Color', ColorTube(count,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
combined_legends{1} = ['Linear shrinkage'];
count = count+1;
% plot for covariance shrinkage
toplot = 2;
param_plot.x_vect       = alpha_vect(alpha_toplot);
P_GP   = PlotTube(Stein_SAA(alpha_toplot, :)', param_plot, ColorTube(toplot,:));
set(P_GP.plot,'facealpha',Trans(toplot));

% add legend
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', 15, 'Location', 'NorthWest')

ylim([250, 500]);
xlim([1e-3, 2e-2]);
ylabel('Stein''s loss','Interpreter','latex','FontSize',FontSize2);
xlabel('$\alpha$','Interpreter','latex','FontSize',FontSize2);
%title('Performance of 3 types of estimator','Interpreter','latex','FontSize',FontSize2);
set(gca,'ygrid','on')
%set(gca,'yscale','log')
set(gca,'xscale','log')

% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%%

beta_toplot  = 1:length(beta_vect);  % QUIC

figure;
hold on; 
combined_legends = [];
count = 3;
% third, plot QUIC
param_plot.x_vect       = beta_vect(beta_toplot);
plot(param_plot.x_vect, average_Stein_QUIC(beta_toplot), LineMarkType{count}, 'Color', ColorTube(count,:), 'MarkerSize', MarkerSize, 'LineWidth', LineWidth);
combined_legends{1} = ['$\ell_1$-regularized ML'];
count = count+1;
% plot tube for QUIC
toplot = 3;
param_plot.x_vect       = beta_vect(beta_toplot);
P_GP   = PlotTube(Stein_QUIC(beta_toplot, :)', param_plot, ColorTube(toplot,:));
set(P_GP.plot,'facealpha',Trans(toplot));

% add legend
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', 15, 'Location', 'NorthWest')

ylim([250, 500]);
ylabel('Stein''s loss','Interpreter','latex','FontSize',FontSize2);
xlabel('$\beta$','Interpreter','latex','FontSize',FontSize2);
%title('Performance of 3 types of estimator','Interpreter','latex','FontSize',FontSize2);
set(gca,'ygrid','on')
%set(gca,'yscale','log')
set(gca,'xscale','log')

% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];