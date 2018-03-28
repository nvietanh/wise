%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasserstein Inverse covariance Shrinkage Estimator
% Viet Anh NGUYEN, Daniel KUHN, Peyman MOHAJERIN ESFAHANI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script plot Figure 1 in the paper.


% Run install.m before running this script
clear all; close all; clc;

%%
rng(0);
S = randn(5);
S = S*S';
[V, D] = eig(S);
d = [0.01, 0.1, 1, 10, 100];
D = diag(d);
S = V*D*V';

rho_vect = 10.^[-1:0.1:2];
eig_tol = 1e-6;
bisection_tol = 1e-6;
for i = 1:length(rho_vect)
    rho = rho_vect(i);
    est = wise( S, rho, eig_tol, bisection_tol);
    eig_value(:, i) = est.eig;
    gamma_value(i) = est.gamma;
    min_gamma_value(i) = est.min_gamma;
    max_gamma_value(i) = est.max_gamma;
end
%%
FontSize = 15;
LineSize = 1.2;
figure;
hold on;
yyaxis left
for i = 1:size(eig_value, 1)
    plot(rho_vect, eig_value(i, :)', '-', 'LineWidth', LineSize);
end
set(gca, 'YScale', 'log')
ylabel('$x^\star$', 'Interpreter', 'latex');

yyaxis right
plot(rho_vect, eig_value(1, :)./eig_value(5, :), '-.o', 'LineWidth', LineSize);
set(gca, 'YScale', 'log')
ylabel('Condition number of $X^\star$', 'Interpreter', 'latex');
set(gca, 'XScale', 'log')
xlabel('$\rho$', 'Interpreter', 'latex');

set(gca,'fontsize',FontSize)

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
FontSize = 15;
LineSize = 1.2;
figure;
hold on;
plot(rho_vect, max_gamma_value, '-.o', 'LineWidth', LineSize, 'Color', [0.8500, 0.3250, 0.0980]);
plot(rho_vect, gamma_value, 'LineWidth', LineSize, 'Color', [0, 0.4470, 0.7410]);
plot(rho_vect, min_gamma_value, 'm-*', 'LineWidth', LineSize);
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('$\rho$', 'Interpreter', 'latex');
ylabel('$\gamma^\star$', 'Interpreter', 'latex');
set(gca,'fontsize',FontSize)
ylim([1e-4 1e3]);

leg_names = {'$\gamma_{\max}$', '$\gamma^\star$', '$\gamma_{\min}$'};
%legend(leg_names, 'Interpreter', 'latex')
legend(leg_names, 'Interpreter', 'latex', 'FontSize', 20)

% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];