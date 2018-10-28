% To plot for Figure 2 and 3, first, load the data
% d50 corresponds to d = 12.5% in the paper
% d200 corresponds to d = 50% in the paper
% d400 corresponds to d = 100% in the paper

% Set plotting parameters

TestCases = 1:100;
% parameters
FontSize  = 24; FontSize2  = 14;
LineWidth = 2;  LineWidth2 = .1;
MarkerSize = 8;

Trans = [0.15, 0.25, 0.35, 0.45, 0.55];



param_plot.Interval = [20,  80];

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
ColorTube(6,:)  = cc(6,:);

LineMarkType = {'-.', ':k', '--', ':', '.', '*'};

y_axis_common = 10.^[0.3 2];
SampleSize_idx_toplot = [1, 2, 3, 4];
%% Plot for Wasserstein
rho_toplot = 5:length(param.rho_vect);
param_plot.x_vect       = param.rho_vect(rho_toplot);

figure;
hold on; 
combined_legends = [];
for toplot = 1:length(SampleSize_idx_toplot)
    plot(param_plot.x_vect, mean(squeeze(AllError.Wasserstein_Robust_Error(SampleSize_idx_toplot(toplot), rho_toplot, TestCases)),2),LineMarkType{toplot},'Color',ColorTube(toplot,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    combined_legends{toplot} = ['$n =$ ' num2str(param.DataSampleSize_vect(SampleSize_idx_toplot(toplot)))];
end


for toplot = 1:length(SampleSize_idx_toplot)
    P_GP   = PlotTube(squeeze(AllError.Wasserstein_Robust_Error(SampleSize_idx_toplot(toplot), rho_toplot, TestCases))' , param_plot, ColorTube(toplot,:));
    set(P_GP.plot,'facealpha',Trans(toplot));
end
xlabel('$\rho$','Interpreter','latex','FontSize',FontSize2);
ylabel('Stein''s Loss','Interpreter','latex','FontSize',FontSize2);
xlim([min(param.rho_vect(rho_toplot)), max(param.rho_vect(rho_toplot))]);
ylim(y_axis_common);
set(gca,'FontSize', FontSize2);
set(gca,'xscale','log')
set(gca,'yscale','log')
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2, 'Location', 'NorthWest')
% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%set(gcf, 'Position', [100, 100, 900, 900])




%% Plot for Linear Shrinkage
alpha_toplot            = 2:(find(param.alpha_vect == 1));
param_plot.x_vect       = param.alpha_vect(alpha_toplot);

figure;
hold on; 
combined_legends = [];
for toplot = 1:length(SampleSize_idx_toplot)
    plot(param_plot.x_vect, mean(squeeze(AllError.Covariance_Shrinkage_Error(SampleSize_idx_toplot(toplot), alpha_toplot, TestCases)),2),'-.','Color',ColorTube(toplot,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    combined_legends{toplot} = ['$n =$ ' num2str(param.DataSampleSize_vect(SampleSize_idx_toplot(toplot)))];
end

for toplot = 1:length(SampleSize_idx_toplot)
    P_GP   = PlotTube(squeeze(AllError.Covariance_Shrinkage_Error(SampleSize_idx_toplot(toplot), alpha_toplot, TestCases))' , param_plot, ColorTube(toplot,:));
    set(P_GP.plot,'facealpha',Trans(toplot));
end
xlabel('$\alpha$','Interpreter','latex','FontSize',FontSize2);
ylabel('Stein''s Loss','Interpreter','latex','FontSize',FontSize2);
xlim([min(param.alpha_vect(alpha_toplot)), max(param.alpha_vect(alpha_toplot))]);
ylim(y_axis_common);
set(gca,'FontSize', FontSize2);
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2, 'Location', 'SouthEast')
set(gca,'xscale','log')
set(gca,'yscale','log')
% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
xticks(10.^[-5 -4 -3 -2 -1 0])
xticklabels({'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^0'})
%% Plot for l1-regularization
beta_toplot = 20:length(param.beta_vect);
param_plot.x_vect = param.beta_vect(beta_toplot);

figure;
ax = axes; hold on; 
ylim(y_axis_common);
ylim manual
set(gca,'yscale','log')
combined_legends = [];

toplot = 1;
for toplot = 1:length(SampleSize_idx_toplot)
    plot(param_plot.x_vect, mean(squeeze(AllError.l1_Error_SDPT3(SampleSize_idx_toplot(toplot), beta_toplot, TestCases)),2),'-.','Color',ColorTube(toplot,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    combined_legends{toplot} = ['$n =$ ' num2str(param.DataSampleSize_vect(SampleSize_idx_toplot(toplot)))];
end

toplot = 1;
for toplot = 1:length(SampleSize_idx_toplot)
    P_GP   = PlotTube(squeeze(AllError.l1_Error_SDPT3(SampleSize_idx_toplot(toplot), beta_toplot, TestCases))' , param_plot, ColorTube(toplot,:));
    set(P_GP.plot,'facealpha',Trans(toplot));
end
xlabel('$\beta$','Interpreter','latex','FontSize',FontSize2);
ylabel('Stein''s Loss','Interpreter','latex','FontSize',FontSize2);
xlim([min(param.beta_vect(beta_toplot)), max(param.beta_vect(beta_toplot))]);
set(gca,'FontSize', FontSize2);
set(gca,'FontSize', FontSize2);
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2, 'Location', 'NorthWest')
set(gca,'xscale','log')
set(gca,'yscale','log')
%ylim manual
% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];



%% The code below plots Figure 3 for the Wasserstein estimator with sparsity structure
% plot Wasserstein structure

k = 3; % change k = 1, 2, 3 for 50%, 75% and 100% of correct information

rho_toplot = 5:length(param.rho_vect);
param_plot.x_vect       = param.rho_vect(rho_toplot);

figure;
hold on; 
combined_legends = [];
for toplot = 1:length(SampleSize_idx_toplot)
    plot(param_plot.x_vect, mean(squeeze(AllError.Wasserstein_Structure_Error(k, SampleSize_idx_toplot(toplot), rho_toplot, TestCases)),2)',LineMarkType{toplot},'Color',ColorTube(toplot,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    combined_legends{toplot} = ['$n =$ ' num2str(param.DataSampleSize_vect(SampleSize_idx_toplot(toplot)))];
end
for toplot = 1:length(SampleSize_idx_toplot)
    P_GP   = PlotTube(squeeze(AllError.Wasserstein_Structure_Error(k, SampleSize_idx_toplot(toplot), rho_toplot, TestCases))' , param_plot, ColorTube(toplot,:));
    set(P_GP.plot,'facealpha',Trans(toplot));
end
xlabel('$\rho$','Interpreter','latex','FontSize',FontSize2);
ylabel('Stein''s Loss','Interpreter','latex','FontSize',FontSize2);
xlim([min(param.rho_vect(rho_toplot)), max(param.rho_vect(rho_toplot))]);
ylim(10.^[0 2]);
set(gca,'FontSize', FontSize2);
%title(['Wasserstein Robust Estimator - ' num2str(max(TestCases)) ' test cases'],'Interpreter','latex','FontSize',FontSize2);
set(gca,'xscale','log')
set(gca,'yscale','log')
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2, 'Location', 'NorthWest')
% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];