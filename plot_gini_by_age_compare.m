function plot_gini_by_age_compare(stats_bench, stats_reform, p, label_bench, label_reform)
%PLOT_GINI_BY_AGE_COMPARE  Plot wealth and disposable income Gini by age.
%
%   plot_gini_by_age_compare(STATS_BENCH, STATS_REFORM, P, LABEL_BENCH, LABEL_REFORM)
%   generates two figures comparing the age-specific Gini coefficients for
%   wealth and disposable income across two scenarios.  The benchmark
%   scenario is shown with a solid line and the tax‑reform scenario with
%   a dashed line.  Each plot has age on the x-axis and the Gini
%   coefficient on the y-axis.  Legends identify the scenarios.

% Extract Gini vectors by age from the stats structures.  These vectors
% should have length p.n_age.
wealth_bench = stats_bench.wealth_gini_by_age;
disp_inc_bench = stats_bench.disp_income_gini_by_age;
wealth_reform = stats_reform.wealth_gini_by_age;
disp_inc_reform = stats_reform.disp_income_gini_by_age;

% X-axis values: ages 1..n_age
ages = 1:p.n_age;

% Plot wealth gini by age
figure; clf; hold on;
plot(ages, wealth_bench, 'b-', 'LineWidth', 2);
plot(ages, wealth_reform, 'r--', 'LineWidth', 2);
hold off;
xlabel('Age');
ylabel('Wealth Gini');
title('Wealth Gini by Age');
legend({label_bench, label_reform}, 'Location','best');
grid on;

% Plot disposable income gini by age
figure; clf; hold on;
plot(ages, disp_inc_bench, 'b-', 'LineWidth', 2);
plot(ages, disp_inc_reform, 'r--', 'LineWidth', 2);
hold off;
xlabel('Age');
ylabel('Disposable income Gini');
title('Disposable Income Gini by Age');
legend({label_bench, label_reform}, 'Location','best');
grid on;

end