function [stats, avgs_by_age_out, pack] = solve_model(p, y_bar_vals, s_next_vals, results_opt)

% SOLVE_MODEL   Steady-state solver for the life-cycle Huggett model.
%
%   [stats, avgs_by_age] = SOLVE_MODEL(p, y_bar_vals, s_next_vals, results_opt) iterates
%   on the aggregate capital stock K and transfer T until the household
%   sector and the firm side are simultaneously consistent.  It then prints
%   macro aggregates plus wealth- and earnings-Gini coefficients.
%
%   INPUTS
%     p            : struct with calibrated parameters
%     y_bar_vals   : n_age x 1 deterministic efficiency profile
%     s_next_vals  : n_age x 1 survival probabilities
%     results_opt  : 1 => generate plots & save .mat, 0 => skip
%
%   OUTPUT
%     stats        : struct with macro and distributional statistics
%     avgs_by_age  : struct containing average consumption, assets, etc. by age
%
% ----------------------------------------------------------------------
% Kaiwen Pi  ·  2025-08-08
% ----------------------------------------------------------------------

% Initialise optional output.  If the caller requests a second output
% argument, this structure will be populated with per-age averages later.
avgs_by_age_out = struct();

%% 0. Demographics
n_z   = p.n_z;
n_age = p.n_age;
R_age = p.R_age;
mu    = calc_age_distribution(p, s_next_vals);  % age weights (sum = 1)
if ~isfield(p, 'G') || isempty(p.G)
    p.G = 0;
end
%% 1. Grids and labour supply
[a_grid, z_grid, a_grid2, Pz, Pz_cdf, Pz1, Pz1_cdf] = ...
    construct_grids_huggett(p);
L = calc_agg_labor(n_z, n_age, z_grid, y_bar_vals, mu, Pz1, Pz);

%% 2. Fixed-point iteration on (K , T_b)
K   = p.K0;    % initial guess for capital
T_b = p.T0;    % initial guess for accidental bequest transfer (T_b)
iter     = 0;
MAX_ITER = 500;
TOL      = 1e-3;    % relative tolerance
wgt      = 0.5;     % damping

while true
    % 2.1 Prices from firm first-order conditions
    w = (1 - p.alpha) * p.TFP * (K / L)^p.alpha;
    r =  p.alpha     * p.TFP * (K / L)^(p.alpha - 1) - p.delta;

    % 2.2 PAYG benefit and implied replacement rate
    ss_ben = p.theta * w * L / sum(mu(R_age:n_age));
    %p.rho  = ss_ben / (w * y_bar_vals(R_age - 1));
    %%% LIFETIME-AVERAGE EARNINGS VERSION
    y_det_lifeavg = mean(y_bar_vals(1:R_age-1));   % lifetime deterministic average
    p.rho= ss_ben / (w * y_det_lifeavg);  % implied replacement rate under lifetime average

    % 2.3 Compute government lump‑sum transfer (Ts).  The government budget
    % collects labour and capital income taxes, net of purchases p.G, and
    % returns the surplus as a per‑capita lump‑sum transfer.  Note that K and
    % L are per‑capita stocks.
    Ts = p.taul * (w * L) + p.tauk * (r * K) - p.G;

    % 2.4 Household optimisation (value-function iteration)
    % The value function takes both the government transfer Ts and the
    % accidental bequest transfer T_b as inputs.
    [v, g] = val_fcn_iter_huggett(a_grid, z_grid, Pz, ...
        y_bar_vals, s_next_vals, Ts, T_b, w, r, p);

    % 2.5 Forward simulation gives new capital and accidental bequest transfer
    % (K', T_b').  Note: Ts is held fixed within this iteration.
    [K_new, T_b_new, avgs_by_age, Phi] = calc_dist_huggett( ...
        a_grid, z_grid, a_grid2, Pz_cdf, Pz1_cdf, ...
        mu, y_bar_vals, s_next_vals, g, w, r, ss_ben, Ts, T_b, p);

    % 2.6 Convergence check
    if abs(K_new - K) <= TOL * K && abs(T_b_new - T_b) <= TOL * T_b
        break;   % converged
    end
    if iter >= MAX_ITER
        warning('solve_model:MaxIter', ...
            'Maximum iterations (%d) reached; exiting.', MAX_ITER);
        break;
    end

    % Relaxed update for stability
    K   = wgt * K_new   + (1 - wgt) * K;
    T_b = wgt * T_b_new + (1 - wgt) * T_b;
    iter = iter + 1;
end

%% 3. Aggregate accounting
Y_hat = p.TFP * (K_new^p.alpha) * (L^(1 - p.alpha));
C_hat = sum(mu .* avgs_by_age.cons);   % aggregate consumption (per capita)
S_hat = K_new;                         % steady-state savings = capital

%% 4. Inequality statistics
% Compute the final government transfer Ts based on the converged values of K and L
Ts_final = p.taul * (w * L) + p.tauk * (r * K_new) - p.G;

% Pass both the accidental bequest transfer (T_b) and the government transfer (Ts)
% into the distribution statistics function.
[stats, Phi_all] = calc_dist_stats(a_grid, a_grid2, z_grid, ...
    y_bar_vals, Phi, mu, w, r, K_new, T_b, Ts_final, L, ss_ben, p);

%% Compute consumption Gini.  Consumption is defined as total net resources
% minus next period's asset choice.  On the coarse asset grid we have
% g(a,z,age) as next‐period assets.  The policy function cons_pol is
% constructed analogously to disposable income but subtracts g to get
% consumption.  We then pass cons_pol into calc_policy_dist_stats to
% obtain the overall consumption Gini.
%
% Consumption for an individual with assets a, productivity state z and
% age j is:
%   c = (1 + (1 - p.tauk)*r) * a + (1 - p.theta - p.taul) * w * e_j(z) + ...
%       Ts_final + T_b + b_j - g(a,z,j)

n_a   = p.n_a;
n_z   = p.n_z;
n_age = p.n_age;

% Precompute the asset part of resources (net interest).  This is a
% column vector of length n_a.
inc_asset_part = (1 + (1 - p.tauk) * r) * a_grid(:);

% Initialise consumption policy function
cons_pol = zeros(n_a, n_z, n_age);

for age = 1:n_age
    % Labour efficiency for this age across z states.  Ensure e_vals is
    % a row vector regardless of z_grid orientation by reshaping the
    % exponential term.
    e_vals = (y_bar_vals(age) * exp(z_grid(:)))';
    % Labour part of resources (net wage): row vector 1-by-n_z
    inc_labour_part = (1 - p.theta - p.taul) * w * e_vals;
    % Social‑security benefit: non‑zero only for retirees
    if age < p.R_age
        b_age = 0.0;
    else
        b_age = ss_ben;
    end
    % Consumption for each a,z: asset part + labour part + transfers + b_age
    % minus next period's asset choice g(a,z,age)
    cons_pol(:,:,age) = repmat(inc_asset_part, 1, n_z) + ...
                       repmat(inc_labour_part, n_a, 1) + ...
                       (Ts_final + T_b + b_age) - g(:,:,age);
end

% Construct a consumption distribution and compute the Gini coefficient.
% We mirror the logic of calc_policy_dist_stats here because that helper
% function is local to calc_dist_stats.m and not visible in this scope.
% First, interpolate the consumption policy onto the fine asset grid (a_grid2).
n_a2  = p.n_a2;
pol2 = zeros(n_a2, n_z, n_age);
for age_idx = 1:n_age
    for iz_idx = 1:n_z
        % linear interpolation from coarse to fine grid
        pol2(:, iz_idx, age_idx) = interp1(a_grid, cons_pol(:, iz_idx, age_idx), a_grid2, 'linear');
    end
end

% Next, construct the density of agents at each grid point on the fine grid.
dist_pdf = zeros(n_a2, n_z, n_age);
for age_idx = 1:n_age
    % first mass point
    dist_pdf(1,:,age_idx) = Phi(1,:,age_idx) * mu(age_idx);
    % differences between cumulative distribution points
    dist_pdf(2:end,:,age_idx) = (Phi(2:end,:,age_idx) - Phi(1:end-1,:,age_idx)) .* mu(age_idx);
end

% Flatten the policy values and corresponding masses into vectors.
pol_vec = pol2(:);
pdf_vec = dist_pdf(:);

% Sort the policy values and reorder the masses accordingly.
[pol_sorted, sort_idx] = sort(pol_vec);
pdf_sorted = pdf_vec(sort_idx);

% Compute the Lorenz curve for consumption.
asset_dist = pol_sorted .* pdf_sorted;
A_tot = sum(asset_dist);
M_tot = sum(pdf_sorted);
if A_tot <= 0 || M_tot <= 0
    cons_gini = NaN;
else
    cum_asset = cumsum(asset_dist) / A_tot;
    cum_pop   = cumsum(pdf_sorted) / M_tot;
    % prepend zero to the cumulative arrays for Lorenz computation
    x_lorenz = [0; cum_pop];
    L_lorenz = [0; cum_asset];
    % Gini coefficient: 1 - sum((L_{i}+L_{i-1})*(x_{i}-x_{i-1}))
    cons_gini = 1 - sum((L_lorenz(2:end) + L_lorenz(1:end-1)) .* diff(x_lorenz));
end

%% Compute ex‑ante expected lifetime utility of a newborn (W_value)
% We evaluate the value function at age 1 (index 1) and zero assets, then take
% expectations over the initial productivity distribution Pz1.  A newborn
% enters with zero assets (a0 = 0).  In our grid, the first element of
% a_grid corresponds to a_min = 0.  We identify the index of zero in the
% asset grid to guard against numerical rounding.  Then we sum over
% productivity states weighted by the initial distribution Pz1.  The
% resulting scalar W_value is stored in the stats struct.
[~, a0_idx] = min(abs(a_grid - 0));
% Value function at a0 for each z at age 1 (first age index)
v_newborn = squeeze(v(a0_idx, :, 1));
% Ensure column vector for dot product
v_newborn = v_newborn(:);
W_value = sum(Pz1(:) .* v_newborn);
stats.W_value = W_value;

%% 5. Console summary
% Print a header showing the scenario name if provided.  This makes
% multiple calls to solve_model easier to distinguish in the console.
if isfield(p, 'scenario_name')
    fprintf('\n====== FINAL STEADY STATE: %s ======\n', p.scenario_name);
else
    fprintf('\n====== FINAL STEADY STATE ======\n');
end
fprintf('Iterations               : %d\n', iter);
fprintf('Capital K                : %.6f\n', K_new);
fprintf('Output  Y                : %.6f\n', Y_hat);
fprintf('Consumption C            : %.6f\n', C_hat);
fprintf('Savings    (=K)          : %.6f\n', S_hat);
fprintf('Labour supply L          : %.6f\n', L);
fprintf('Wage w                   : %.6f\n', w);
fprintf('Interest rate r          : %.6f\n', r);
fprintf('Social-security benefit  : %.6f\n', ss_ben);
fprintf('Replacement rate rho     : %.6f\n', p.rho);
fprintf('K / Y                    : %.6f\n', K_new / Y_hat);
fprintf('C / Y                    : %.6f\n', C_hat / Y_hat);
fprintf('Tax rates (tau_k, tau_l) : (%.4f , %.4f)\n', p.tauk, p.taul);
fprintf('\n-- Inequality --\n');
fprintf('Wealth  Gini             : %.4f\n', stats.wealth_gini);
fprintf('Earnings Gini            : %.4f\n', stats.earnings_gini);
fprintf('Disposable income Gini   : %.4f\n', stats.disp_income_gini);
fprintf('Consumption Gini         : %.4f\n', cons_gini);
        % Additional distributional statistics.
        % Print wealth shares by percentile groups.  These statistics are computed
        % in calc_dist_stats.m and capture the cumulative share of total wealth
        % held by different segments of the population.  For example, the
        % "Top 1%% share of wealth" reports the fraction of aggregate wealth owned
        % by the richest 1%% of households.  The bottom shares are the cumulative
        % share up to the indicated percentile.
        fprintf('-- Wealth share by percentile --\n');
        fprintf('Top 1%% share of wealth              : %.4f\n', stats.wealth_top_01);
        fprintf('Top 5%% share of wealth              : %.4f\n', stats.wealth_top_05);
        fprintf('Top 20%% share of wealth             : %.4f\n', stats.wealth_top_20);
        fprintf('Top 40%% share of wealth             : %.4f\n', stats.wealth_top_40);
        fprintf('Top 60%% share of wealth             : %.4f\n', stats.wealth_top_60);
        fprintf('Top 80%% share of wealth             : %.4f\n', stats.wealth_top_80);
        fprintf('Bottom 20%% share of wealth          : %.4f\n', stats.wealth_bottom_20);
        fprintf('Bottom 40%% share of wealth          : %.4f\n', stats.wealth_bottom_40);
        % Print earnings (pre-tax labour income) shares by percentile.  These
        % measures show the cumulative share of total labour income earned by
        % different income strata.
        fprintf('-- Earnings share by percentile --\n');
        fprintf('Top 1%% share of earnings            : %.4f\n', stats.earnings_top_01);
        fprintf('Top 5%% share of earnings            : %.4f\n', stats.earnings_top_05);
        fprintf('Top 20%% share of earnings           : %.4f\n', stats.earnings_top_20);
        fprintf('Top 40%% share of earnings           : %.4f\n', stats.earnings_top_40);
        fprintf('Top 60%% share of earnings           : %.4f\n', stats.earnings_top_60);
        fprintf('Top 80%% share of earnings           : %.4f\n', stats.earnings_top_80);
        fprintf('Bottom 20%% share of earnings        : %.4f\n', stats.earnings_bottom_20);
        fprintf('Bottom 40%% share of earnings        : %.4f\n', stats.earnings_bottom_40);
        % Print disposable income shares by percentile.  Disposable income
        % includes net-of-tax labour and capital income plus transfers and
        % benefits.  These shares illustrate how after-tax resources are
        % distributed across households.
        fprintf('-- Disposable income share by percentile --\n');
        fprintf('Top 1%% share of disposable income   : %.4f\n', stats.disp_income_top_01);
        fprintf('Top 5%% share of disposable income   : %.4f\n', stats.disp_income_top_05);
        fprintf('Top 20%% share of disposable income  : %.4f\n', stats.disp_income_top_20);
        fprintf('Top 40%% share of disposable income  : %.4f\n', stats.disp_income_top_40);
        fprintf('Top 60%% share of disposable income  : %.4f\n', stats.disp_income_top_60);
        fprintf('Top 80%% share of disposable income  : %.4f\n', stats.disp_income_top_80);
        fprintf('Bottom 20%% share of disp. income    : %.4f\n', stats.disp_income_bottom_20);
        fprintf('Bottom 40%% share of disp. income    : %.4f\n', stats.disp_income_bottom_40);
        % Print Gini coefficients explicitly for clarity.  These are measures of
        % inequality for wealth, earnings and disposable income, respectively.

%% 5.1 Government transfer Ts (per‑capita; excludes bequests)
% Budget:  G + Ts = τ_l·w·L + τ_k·r·K.  We report the per‑capita Ts.
% Total population. If you simulate NN per age, N = NN * n_age.
if isfield(p,'NN') && isfield(p,'n_age')
    N_pop = p.NN * p.n_age;
else
    N_pop = 1;  % if everything is already normalized per capita
end

if ~isfield(p,'G'),    p.G    = 0; end    % per-capita government purchases

% Per-capita objects (your code keeps K_new, C_hat per capita)
K_pc = K_new;
G_pc = p.G;


% Per‑capita Ts from the budget constraint
Ts_pc =  p.taul * (w * L) + p.tauk * (r * K_pc) - G_pc;

% Aggregate TOTAL Ts and its share (not printed by default)
Ts_total = N_pop * Ts_pc;
Y_total  = Y_hat * N_pop;

% Print the per‑capita transfer.  This is the Ts that enters the household budget
% constraint.
fprintf('Government Transfer Ts (per capita) : %.6f\n', Ts_pc);
% Uncomment the following lines to print aggregate Ts and its share of output.
% fprintf('Government Transfer Ts (TOTAL) : %.6f\n', Ts_total);
% fprintf('Ts / Y (aggregate)             : %.6f\n', Ts_total / max(1e-12, Y_total));
if Ts_pc < 0
    fprintf('Note: Ts < 0 (head tax implied).\n');
end

%% 6. Optional plots and save
if results_opt ~= 0
    fig_id = 0;
    fig_id = plot_distribution(fig_id, a_grid2, Phi, Phi_all, p);
    fig_id = plot_value_and_policy_functions(fig_id, a_grid, v, g, p);
    fig_id = plot_vars_by_age(fig_id, avgs_by_age, p);
    save('results.mat', 'stats', 'avgs_by_age', 'Phi');
end

% If the caller requested a second output, populate it with the
% age-specific averages computed inside the fixed-point iteration.  When
% solve_model is invoked with only one output, this assignment is
% harmless because MATLAB will simply ignore the extra return values.
if nargout >= 2
    avgs_by_age_out = avgs_by_age;
end
% --- Optional 3rd output: objects needed for welfare comparison ---
if nargout >= 3
    pack = struct();
    pack.v       = v;         % value function on coarse asset grid (n_a x n_z x n_age)
    pack.Phi     = Phi;       % CDF on fine asset grid by age (n_a2 x n_z x n_age)
    pack.mu      = mu;        % age weights
    pack.a_grid  = a_grid;    % coarse asset grid
    pack.a_grid2 = a_grid2;   % fine asset grid
    pack.z_grid  = z_grid;    % productivity grid
end
% --- Optional 3rd output: objects needed for welfare comparison ---
if nargout >= 3
    pack = struct();
    pack.v       = v;         % value function on coarse asset grid (n_a x n_z x n_age)
    pack.Phi     = Phi;       % CDF on fine asset grid by age (n_a2 x n_z x n_age)
    pack.mu      = mu;        % age weights
    pack.a_grid  = a_grid;    % coarse asset grid
    pack.a_grid2 = a_grid2;   % fine asset grid
    pack.z_grid  = z_grid;    % productivity grid
end

end
