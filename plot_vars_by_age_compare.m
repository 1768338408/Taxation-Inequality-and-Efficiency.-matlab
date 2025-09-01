function [fn] = plot_vars_by_age_compare(fn, avgs1, avgs2, p, label1, label2)
%PLOT_VARS_BY_AGE_COMPARE  Compare life-cycle profiles under two scenarios.
%
%   This helper function plots average consumption, capital, next-period
%   capital, earnings and disposable income by age for two sets of results.
%   It overlays the life-cycle profiles in the same figure for easy
%   comparison.  The user must supply descriptive labels for the legend
%   entries corresponding to each scenario.
%
%   INPUTS
%     fn      - starting figure number (incremented internally)
%     avgs1   - struct with fields .cons, .K, .Kn, .earn, .inc for scenario 1
%     avgs2   - struct with the same fields for scenario 2
%     p       - parameter struct; only p.n_age is used here
%     label1  - legend label for scenario 1 (e.g. 'Tax Reform')
%     label2  - legend label for scenario 2 (e.g. 'Benchmark')
%
%   OUTPUT
%     fn      - final figure number after plotting

% number of age cohorts
n_age = p.n_age;
ages  = 1:n_age;

% -----------------------------------------------------------------------
% Compute the stationary age weights (mu) from the input vectors file.
%
% These weights determine the fraction of the population in each age
% cohort.  They are computed from the survival probabilities stored in
% ``input_vectors.txt``.  If the file cannot be read, we fall back to
% uniform weights.
try
    fid_mu = fopen('input_vectors.txt','r');
    if fid_mu >= 0
        data_mu = textscan(fid_mu, '%d %f %f');
        fclose(fid_mu);
        s_next_vals = data_mu{3};
        mu = ones(n_age,1);
        for idx = 2:n_age
            mu(idx) = mu(idx-1) * s_next_vals(idx-1);
        end
        mu = mu / sum(mu);
    else
        mu = ones(n_age,1) / n_age;
    end
catch
    mu = ones(n_age,1) / n_age;
end

% List of variable names and y-axis labels.  The order here determines
% the order in which the figures are generated.
vars   = {'cons','K','Kn','earn','inc'};
ylabs  = {'Consumption','Capital','Next-period capital','Earnings','Disposable income'};

% Loop over each variable to produce an overlay plot.
for i = 1:numel(vars)
    var_name = vars{i};
    ylab     = ylabs{i};

    fn = fn + 1;
    figure(fn); clf;
    hold on;
    % Retrieve per-age series
    data1 = avgs1.(var_name);
    data2 = avgs2.(var_name);
    current_ylab = ylab;
    % For the four economic variables (consumption, capital, earnings and
    % disposable income), aggregate across the entire population rather
    % than displaying per-capita averages.  We multiply the average by
    % ``mu`` (the population share at each age) and by the total number
    % of simulated individuals (``p.NN``) across all age groups.  This
    % yields a total value for each age cohort.  ``Kn`` (next-period
    % capital) remains unscaled since it is a transitional variable.
    if any(strcmp(var_name, {'cons','K','earn','inc'}))
        % Determine total population represented in the simulation
        if isfield(p, 'NN')
            total_pop = p.NN * n_age;
        else
            total_pop = n_age;
        end
        data1 = data1 .* mu * total_pop;
        data2 = data2 .* mu * total_pop;
        % Choose an appropriate scale to avoid excessively large numbers.
        max_val = max([max(data1), max(data2)]);
        if max_val >= 1e4
            scale_factor = 1e4;
            unit_suffix  = '10k RMB';
        elseif max_val >= 1e3
            scale_factor = 1e3;
            unit_suffix  = 'k RMB';
        else
            scale_factor = 1;
            unit_suffix  = 'RMB';
        end
        data1 = data1 / scale_factor;
        data2 = data2 / scale_factor;
        current_ylab = sprintf('%s (%s)', ylab, unit_suffix);
    end
    % Plot scenarios with specified line styles: scenario 1 (tax reform)
    % is dashed red; scenario 2 (benchmark) is solid blue.
    plot(ages, data1, 'r--', 'LineWidth', 2);
    plot(ages, data2, 'b-',  'LineWidth', 2);
    hold off;
    xlabel('Age');
    ylabel(current_ylab);
    legend({label1, label2}, 'Location', 'best');
    title(sprintf('%s by age', ylab));
end

end