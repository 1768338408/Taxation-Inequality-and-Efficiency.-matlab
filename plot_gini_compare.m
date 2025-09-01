function plot_gini_compare(stats_bench, stats_reform, label_bench, label_reform)
%PLOT_GINI_COMPARE Plot wealth and disposable income Gini coefficients.
%
%   plot_gini_compare(STATS_BENCH, STATS_REFORM, LABEL_BENCH, LABEL_REFORM)
%   produces a figure comparing the aggregated wealth Gini and disposable
%   income Gini coefficients for two scenarios.  The x-axis lists the two
%   variables ("Wealth Gini" and "Disposable income Gini"), and the y-axis
%   shows the Gini coefficient (unitless).  The first scenario is drawn
%   using a solid line and the second using a dashed line.  The caller
%   should supply descriptive labels for the legend entries.

% extract gini values
gini_bench = [stats_bench.wealth_gini, stats_bench.disp_income_gini];
gini_reform = [stats_reform.wealth_gini, stats_reform.disp_income_gini];

% x positions for the two variables
x_vals = 1:2;

figure; clf; hold on;
% benchmark: solid line
plot(x_vals, gini_bench, 'b-', 'LineWidth', 2);
% tax reform: dashed line
plot(x_vals, gini_reform, 'r--', 'LineWidth', 2);

% axes labels and formatting
set(gca, 'XTick', x_vals, 'XTickLabel', {'Wealth Gini','Disposable income Gini'});
xlabel('Gini measure');
ylabel('Gini coefficient (unitless)');
title('Aggregated Gini Coefficients by Scenario');
legend({label_bench, label_reform}, 'Location','best');
grid on;
hold off;

end