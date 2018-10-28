
TestCases = 1:100;
LossType = 'Log-likelihood';
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

%LineMarkType = {'-.', ':k', '--', ':', '.', '*'};
LineMarkType = {'-', ':k', '--', '-o', '.', '*'};

%Leg = legend('Wass');
%set(Leg,'Interpreter','latex','FontSize', FontSize2)

%P_SAA = PlotTube(ones(size(param_plot.rho_vect))*SAError(:, SampleSize_idx), param_plot, ColorTube1);
%set(P_SAA.plot,'facealpha',Trans-0.1);

y_axis_common = 10.^[0.3 2];
SampleSize_idx_toplot = [1, 2, 3, 4];




%% Plot learning curve for Wasserstein in another definition of rho
% Here, the optimal rho is the average of optimal rho
files = {'LearningCurve50.mat', 'LearningCurve200.mat', 'LearningCurve400.mat'};


figure;
hold on;
for count = 1:length(files)
    load(files{count}, 'AllError', 'param');
    optimal_rho = [];
    temp = real(AllError.Wasserstein_Robust_Error);
    param_plot.x_vect = param.DataSampleSize_vect;
    for i = 1:length(param_plot.x_vect)
        [val, idx] = min(squeeze(temp(i, :, :)));
        optimal_rho(i) = mean(param.rho_vect(idx));
    end
    plot(param_plot.x_vect, optimal_rho, LineMarkType{count},'Color',ColorTube(count,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);

    % run regression to find the slope of the curve
    display('Regression coefficient:');
    regress(log(optimal_rho)', log(param_plot.x_vect)')
end
combined_legends = {'$d=12.5\%$', '$d=50\%$', '$d=100\%$'};
xlabel('$n$','Interpreter','latex','FontSize',FontSize2);
ylabel('$\rho^*$','Interpreter','latex','FontSize',FontSize2);
xlim([10, 10^6]);
ylim([0, 2]);

xticks(10.^[1:1:6])
xticklabels({'10^1','10^2','10^3','10^4','10^5','10^6'})

set(gca,'FontSize', FontSize2);
%title(['Wasserstein Robust Estimator - ' num2str(max(TestCases)) ' test cases'],'Interpreter','latex','FontSize',FontSize2);
Leg = legend(combined_legends);
set(Leg,'Interpreter','latex','FontSize', FontSize2+5, 'Location', 'NorthEast')

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