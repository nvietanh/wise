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

%Leg = legend('Wass');
%set(Leg,'Interpreter','latex','FontSize', FontSize2)

%P_SAA = PlotTube(ones(size(param_plot.rho_vect))*SAError(:, SampleSize_idx), param_plot, ColorTube1);
%set(P_SAA.plot,'facealpha',Trans-0.1);

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
%set(gcf, 'Position', [100, 100, 900, 900])


%% plot Wasserstein structure

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

%% Plot for Linear Shrinkage
alpha_toplot            = 1:(find(param.alpha_vect == 1));
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
%title(['Shrinkage Estimator - ' num2str(max(TestCases)) ' test cases'],'Interpreter','latex','FontSize',FontSize2);
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
%title(['$l_1$-Regularized Estimator - ' num2str(max(TestCases)) ' test cases'],'Interpreter','latex','FontSize',FontSize2);

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


%% Plot for Ledoit Wolf shrinkage
alpha_toplot            = 1:length(param.alpha_vect);
param_plot.x_vect       = param.alpha_vect(alpha_toplot);

figure;
hold on; 
combined_legends = [];
for toplot = 1:length(SampleSize_idx_toplot)
    plot(param_plot.x_vect, mean(squeeze(AllError.LedoitWolf_Shrinkage_Error(SampleSize_idx_toplot(toplot), alpha_toplot, TestCases)),2),'-.','Color',ColorTube(toplot,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    combined_legends{toplot} = ['$n =$ ' num2str(param.DataSampleSize_vect(SampleSize_idx_toplot(toplot)))];
end

for toplot = 1:length(SampleSize_idx_toplot)
    P_GP   = PlotTube(squeeze(AllError.LedoitWolf_Shrinkage_Error(SampleSize_idx_toplot(toplot), alpha_toplot, TestCases))' , param_plot, ColorTube(toplot,:));
    set(P_GP.plot,'facealpha',Trans(toplot));
end
xlabel('$\alpha$','Interpreter','latex','FontSize',FontSize2);
ylabel(['Error'],'Interpreter','latex','FontSize',FontSize2);
xlim([min(param.alpha_vect(alpha_toplot)), max(param.alpha_vect(alpha_toplot))]);
ylim(y_axis_common);
set(gca,'FontSize', FontSize2);
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2, 'Location', 'SouthEast')
set(gca,'xscale','log')
set(gca,'yscale','log')
%title(['Shrinkage Estimator - ' num2str(max(TestCases)) ' test cases'],'Interpreter','latex','FontSize',FontSize2);
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

%% Plot all 4 estimators in 1 plot
LossType = 'Log-likelihood';

% only one value for the samplesize!!!
SampleSize_idx_toplot = 3;
rho_toplot = 50:length(param.rho_vect); % Wasserstein
alpha_toplot = 1:(find(param.alpha_vect == 1)); % shrinkage
beta_toplot = 30:(length(param.beta_vect)); % l-1


% set toplot = 1; 
toplot = 1;

figure;
hold on; 
combined_legends = [];

count = 1;
% first, plot Wasserstein
param_plot.x_vect       = param.rho_vect(rho_toplot);
plot(param_plot.x_vect, mean(squeeze(AllError.Wasserstein_Robust_Error(SampleSize_idx_toplot, rho_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
combined_legends{count} = ['Wasserstein Robust'];
count = count+1;

% second, plot Covariance Shrinkage
param_plot.x_vect       = param.alpha_vect(alpha_toplot);
plot(param_plot.x_vect, mean(squeeze(AllError.Covariance_Shrinkage_Error(SampleSize_idx_toplot, alpha_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
combined_legends{count} = ['Covariance Shrinkage'];
count = count+1;

% finally, plot l1-estimator
param_plot.x_vect       = param.beta_vect(beta_toplot);
plot(param_plot.x_vect, mean(squeeze(AllError.DAspremont_Error(SampleSize_idx_toplot, beta_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
combined_legends{count} = ['$l_1$'];
count = count+1;

% plot wasserstein structure
for k = 1:length(param.zero_mat)
    param_plot.x_vect       = param.rho_vect(rho_toplot);
    plot(param_plot.x_vect, mean(squeeze(AllError.Wasserstein_Structure_Error(k, SampleSize_idx_toplot, rho_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    combined_legends{count} = ['Wasserstein Structure ' num2str(pct_info(k)*100)];
    count = count+1;
end


% now, start to plot the tube
count = 1;

param_plot.x_vect       = param.rho_vect(rho_toplot);
P_GP   = PlotTube(squeeze(AllError.Wasserstein_Robust_Error(SampleSize_idx_toplot(toplot), rho_toplot, TestCases))' , param_plot, ColorTube(count,:));
set(P_GP.plot,'facealpha',Trans(count));
count = count+1;

param_plot.x_vect       = param.alpha_vect(alpha_toplot);
P_GP   = PlotTube(squeeze(AllError.Covariance_Shrinkage_Error(SampleSize_idx_toplot(toplot), alpha_toplot, TestCases))' , param_plot, ColorTube(count,:));
set(P_GP.plot,'facealpha',Trans(count));
count = count+1;


param_plot.x_vect       = param.beta_vect(beta_toplot);
P_GP   = PlotTube(squeeze(AllError.DAspremont_Error(SampleSize_idx_toplot(toplot), beta_toplot, TestCases))' , param_plot, ColorTube(count,:));
set(P_GP.plot,'facealpha',Trans(count));

ylabel(['Error - Measured by ' LossType ' Loss'],'Interpreter','latex','FontSize',FontSize2);
title(['Performance of 4 types of estimator - ' num2str(max(TestCases)) ' test cases, p = 20, n = ' num2str(param.DataSampleSize_vect(SampleSize_idx_toplot(toplot)))],'Interpreter','latex','FontSize',FontSize2);
set(gca,'ygrid','on')

set(gca,'xscale','log')
set(gca,'yscale','log')
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2)



%% Plot learning curve for Wasserstein
% Here, the optimal rho is the average of the optimal rho for each test
% case
files = {'sparse_d400_LearningCurve2.mat', 'sparse_d200_LearningCurve2.mat', 'sparse_d100_LearningCurve2.mat', 'sparse_d50_LearningCurve2.mat'};


figure;
hold on;
for count = 1:4
    load(['resultdata/' files{count}], 'AllError', 'param');
    temp = real(AllError.Wasserstein_Robust_Error);
    param_plot.x_vect = param.DataSampleSize_vect;
    for i = 1:length(param_plot.x_vect)
        [~, idx] = min(squeeze(temp(i, :, :)));
        y(i) = mean(param.rho_vect(idx));
    end
    plot(param_plot.x_vect, y, LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    haha{count} = temp;
    hehe{count} = y;
end
combined_legends = {'$d=400$', '$d=200$', '$d=100$', '$d=50$'};
xlabel('$n$','Interpreter','latex','FontSize',FontSize2);
ylabel('$\rho^*$','Interpreter','latex','FontSize',FontSize2);
xlim([min(param_plot.x_vect), max(param_plot.x_vect)]);
ylim([0, 20]);
set(gca,'FontSize', FontSize2);
%title(['Wasserstein Robust Estimator - ' num2str(max(TestCases)) ' test cases'],'Interpreter','latex','FontSize',FontSize2);
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2, 'Location', 'NorthEast')

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

%%
% Plot all 3 estimators in 1 plot
LossType = 'Log-likelihood';

% only one value for the samplesize!!!
SampleSize_idx_toplot = [1, 2, 3, 4];
rho_toplot = 1:length(param.rho_vect); % Wasserstein
alpha_toplot = 1:(find(param.alpha_vect == 1)); % shrinkage
beta_toplot = 1:(length(param.beta_vect)); % l-1


combined_legends = {'Wasserstein Robust', 'Covariance Shrinkage', '$l_1$'};
count = 4;
for i = 1:length(param.zero_mat)
        combined_legends{count} = ['Wasserstein Structure ' num2str(pct_info(i)*100)];
        count = count+1;
end
figure;
hold on; 

for i = 1:length(SampleSize_idx_toplot)
    subplot(2, 2, i);
    hold on
    count = 1;
    % first, plot Wasserstein
    param_plot.x_vect       = param.rho_vect(rho_toplot);
    plot(param_plot.x_vect, mean(squeeze(AllError.Wasserstein_Robust_Error(SampleSize_idx_toplot(i), rho_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    count = count+1;

    % second, plot Covariance Shrinkage
    param_plot.x_vect       = param.alpha_vect(alpha_toplot);
    plot(param_plot.x_vect, mean(squeeze(AllError.Covariance_Shrinkage_Error(SampleSize_idx_toplot(i), alpha_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    count = count+1;

    % finally, plot l1-estimator
    param_plot.x_vect       = param.beta_vect(beta_toplot);
    plot(param_plot.x_vect, mean(squeeze(AllError.DAspremont_Error(SampleSize_idx_toplot(i), beta_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);

    count = count+1;

    % plot wasserstein structure
    for k = 1:length(param.zero_mat)
        param_plot.x_vect       = param.rho_vect(rho_toplot);
        plot(param_plot.x_vect, mean(squeeze(AllError.Wasserstein_Structure_Error(k, SampleSize_idx_toplot(i), rho_toplot, TestCases)),2),LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
        count = count + 1;
    end


    % now, start to plot the tube
    count = 1;

    param_plot.x_vect       = param.rho_vect(rho_toplot);
    P_GP   = PlotTube(squeeze(AllError.Wasserstein_Robust_Error(SampleSize_idx_toplot(i), rho_toplot, TestCases))' , param_plot, ColorTube(count,:));
    set(P_GP.plot,'facealpha',Trans(count));
    count = count+1;

    param_plot.x_vect       = param.alpha_vect(alpha_toplot);
    P_GP   = PlotTube(squeeze(AllError.Covariance_Shrinkage_Error(SampleSize_idx_toplot(i), alpha_toplot, TestCases))' , param_plot, ColorTube(count,:));
    set(P_GP.plot,'facealpha',Trans(count));
    count = count+1;


    param_plot.x_vect       = param.beta_vect(beta_toplot);
    P_GP   = PlotTube(squeeze(AllError.DAspremont_Error(SampleSize_idx_toplot(i), beta_toplot, TestCases))' , param_plot, ColorTube(count,:));
    set(P_GP.plot,'facealpha',Trans(count));

    ylabel(['Error'],'Interpreter','latex','FontSize',FontSize2);
    title(['n = ' num2str(param.DataSampleSize_vect(SampleSize_idx_toplot(i)))],'Interpreter','latex','FontSize',FontSize2);
    set(gca,'ygrid','on')

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1, 100]);
    xlim([10^-2, 10]);
end
Leg = legend(combined_legends);
Leg.Orientation = 'horizontal';
set(Leg,'Interpreter','latex','FontSize', FontSize2)