function main

close all

% load parameters
p = load_parameters('input_params.txt');

% load other inputs
[ y_bar_vals, s_next_vals ] = load_other_inputs( 'input_vectors.txt', p.n_age );

if (p.run_opt == 1)
    % Compute both the tax reform and benchmark steady states. The first
    % scenario uses the tax rates read from the input file. The second
    % scenario uses benchmark tax rates (tau_k=0.165, tau_l=0.03). We
    % disable internal plotting inside solve_model and instead overlay
    % the life-cycle profiles of the key variables using a dedicated
    % comparison plotting routine.

    % Reform scenario (uses p.tauk and p.taul from input)
    p_reform = p;
    p_reform.scenario_name = 'Tax Reform';
    results_opt = 0;
    [stats_reform, avgs_reform] = solve_model(p_reform, y_bar_vals, s_next_vals, results_opt);

    % Benchmark scenario
    p_benchmark = p;
    p_benchmark.tauk = 0.165;
    p_benchmark.taul = 0.03;
    p_benchmark.scenario_name = 'Benchmark';
    [stats_benchmark, avgs_benchmark] = solve_model(p_benchmark, y_bar_vals, s_next_vals, results_opt);

    % Overlay key life-cycle profiles
    fig_id = 0;
    fig_id = plot_vars_by_age_compare(fig_id, avgs_reform, avgs_benchmark, p, 'Tax Reform', 'Benchmark');

    % Summary print
    fprintf('\n===== Comparison Summary =====\n');
    fprintf('Tax Reform scenario:    tau_k = %.3f, tau_l = %.3f\n', p_reform.tauk, p_reform.taul);
    fprintf('Benchmark scenario:     tau_k = %.3f, tau_l = %.3f\n', p_benchmark.tauk, p_benchmark.taul);
    fprintf('Capital K - reform vs. benchmark:     %.6f  vs  %.6f\n', stats_reform.K, stats_benchmark.K);
    fprintf('Output  Y - reform vs. benchmark:     %.6f  vs  %.6f\n', stats_reform.Y, stats_benchmark.Y);
    fprintf('Wealth Gini - reform vs. benchmark:    %.4f  vs  %.4f\n', stats_reform.wealth_gini, stats_benchmark.wealth_gini);
    fprintf('Earnings Gini - reform vs. benchmark:  %.4f  vs  %.4f\n', stats_reform.earnings_gini, stats_benchmark.earnings_gini);
    fprintf('Disposable income Gini - reform vs. benchmark: %.4f  vs  %.4f\n', stats_reform.disp_income_gini, stats_benchmark.disp_income_gini);

    plot_gini_compare(stats_benchmark, stats_reform, 'Benchmark', 'Tax Reform');
    plot_gini_by_age_compare(stats_benchmark, stats_reform, p, 'Benchmark', 'Tax Reform');

    % ---------- Overall CEV for newborn (ex-ante) ----------
    W_R = stats_reform.W_value;
    W_B = stats_benchmark.W_value;
    sigma = p.sigma;
    if abs(sigma - 1) < 1e-8
        CEV_val = exp(W_R - W_B) - 1;
    else
        CEV_val = (W_R / W_B)^(1/(1 - sigma)) - 1;
    end
    fprintf('Consumption-equivalent variation (CEV, reform vs. benchmark) : %.6f\n', CEV_val);

    % ===== Pct. of Agents Better Off (weighted by benchmark cross-section) =====
try
    % 1) Grids and Markov objects
    [a_grid, z_grid, a_grid2, Pz, Pz_cdf, ~, Pz1_cdf] = construct_grids_huggett(p_reform);

    % 2) Set rho consistent with printed replacement (from ss_ben)
    R_age = p.R_age;
    y_det_lifeavg = mean(y_bar_vals(1:R_age-1));
    pR = p_reform;    pR.rho = stats_reform.ss_ben    / (stats_reform.w    * y_det_lifeavg);
    pB = p_benchmark; pB.rho = stats_benchmark.ss_ben / (stats_benchmark.w * y_det_lifeavg);

    % 3) One-shot value function iteration at fixed prices/transfers
    [vR, gR] = val_fcn_iter_huggett( ...
        a_grid, z_grid, Pz, y_bar_vals, s_next_vals, ...
        stats_reform.Ts, stats_reform.T_b, stats_reform.w, stats_reform.r, pR);

    [vB, gB] = val_fcn_iter_huggett( ...
        a_grid, z_grid, Pz, y_bar_vals, s_next_vals, ...
        stats_benchmark.Ts, stats_benchmark.T_b, stats_benchmark.w, stats_benchmark.r, pB);

    % 4) Benchmark cross-sectional distribution (CDF) for weighting
    mu = calc_age_distribution(p, s_next_vals);   % age shares
    [~, ~, ~, Phi_B] = calc_dist_huggett( ...
        a_grid, z_grid, a_grid2, Pz_cdf, Pz1_cdf, ...
        mu, y_bar_vals, s_next_vals, gB, ...
        stats_benchmark.w, stats_benchmark.r, stats_benchmark.ss_ben, ...
        stats_benchmark.Ts, stats_benchmark.T_b, pB);

    % 5) Interpolate values to fine asset grid
    n_a2 = length(a_grid2); n_z = length(z_grid); n_age = p.n_age;
    vR_fine = zeros(n_a2, n_z, n_age);
    vB_fine = zeros(n_a2, n_z, n_age);
    for age = 1:n_age
        for iz = 1:n_z
            vR_fine(:,iz,age) = interp1(a_grid, vR(:,iz,age), a_grid2, 'linear', 'extrap');
            vB_fine(:,iz,age) = interp1(a_grid, vB(:,iz,age), a_grid2, 'linear', 'extrap');
        end
    end

    % 6) CDF -> PDF and apply age weights mu
    dist_pdf = zeros(n_a2, n_z, n_age);
    for age = 1:n_age
        dist_pdf(1,:,age)     = Phi_B(1,:,age) * mu(age);
        dist_pdf(2:end,:,age) = (Phi_B(2:end,:,age) - Phi_B(1:end-1,:,age)) * mu(age);
    end

    % 7) Winners mass
    winners_mask = (vR_fine - vB_fine) >= 0;
    mass_total   = sum(dist_pdf(:));
    mass_winners = sum(dist_pdf(winners_mask));
    pct_winners  = 100 * mass_winners / max(1e-16, mass_total);
    fprintf('Pct. of Agents Better Off (reform vs. benchmark): %.1f\n', pct_winners);

    % ---------- CEV at birth by z ----------
    try
        [~, a0_idx_birth] = min(abs(a_grid - 0));
        vR_birth_z  = squeeze(vR(a0_idx_birth, :, 1));
        vB_birth_z  = squeeze(vB(a0_idx_birth, :, 1));
        sigma_z = p.sigma;
        if abs(sigma_z - 1) < 1e-8
            cev_birth_by_z = exp(vR_birth_z - vB_birth_z) - 1;
        else
            cev_birth_by_z = (vR_birth_z ./ vB_birth_z).^(1/(1 - sigma_z)) - 1;
        end
        fprintf('\n=== CEV at birth by z ===\n');
        fprintf('  (a0 = 0, age j = 1; sigma = %.4f)\n', sigma_z);
        for iz = 1:length(z_grid)
            fprintf('  z[%2d] = %+0.6f : CEV_birth = %+0.6f\n', iz, z_grid(iz), cev_birth_by_z(iz));
        end

        % ---------- Weighted CEV by terciles (initial z PDF) ----------
        try
            [~, ~, ~, ~, ~, Pz1, ~] = construct_grids_huggett(p_reform);
            weights_z = Pz1(:);
            cev_z     = cev_birth_by_z(:);
            n_z_total = length(z_grid);
            group_size = floor(n_z_total/3);
            idx_bottom = 1:group_size;
            idx_middle = (group_size+1):(2*group_size);
            idx_top    = (2*group_size+1):n_z_total;

            w_bottom = weights_z(idx_bottom);  w_middle = weights_z(idx_middle);  w_top = weights_z(idx_top);
            c_bottom = cev_z(idx_bottom);      c_middle = cev_z(idx_middle);      c_top = cev_z(idx_top);

            wcev_bottom = sum(w_bottom .* c_bottom) / max(1e-16, sum(w_bottom));
            wcev_middle = sum(w_middle .* c_middle) / max(1e-16, sum(w_middle));
            wcev_top    = sum(w_top    .* c_top)    / max(1e-16, sum(w_top));

            fprintf('\n=== CEV at birth by productivity group (weighted) ===\n');
            fprintf('  Bottom 1/3 weighted CEV_birth: %+0.6f\n', wcev_bottom);
            fprintf('  Middle 1/3 weighted CEV_birth: %+0.6f\n', wcev_middle);
            fprintf('  Top 1/3 weighted CEV_birth: %+0.6f\n', wcev_top);
        catch group_err
            warning('WeightedCEV:Failed','Failed to compute weighted CEV by productivity group: %s', group_err.message);
        end
    catch cevErr
        warning('CEVbyZ:Failed','Failed to compute CEV by z: %s', cevErr.message);
    end
catch ME
    warning('PctBetterOff:Failed','Failed to compute Pct. Better Off: %s', ME.message);
end


elseif (p.run_opt == 2)
    % calibrate model
    calibrate_model(p, y_bar_vals, s_next_vals);
else
    error('invalid value for run_opt')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ p ] = load_parameters( filename )

p = {};

% fopen returns fid >= 0 if successful, -1 if not
fid = fopen(filename, 'r');

% read data, assume file is valid
data_cell = textscan(fid, '%f%*[^\n]');
data_array = data_cell{1};

i = 0;
i=i+1; p.n_a        = data_array(i);
i=i+1; p.n_z        = data_array(i);
i=i+1; p.n_age      = data_array(i);
i=i+1; p.R_age      = data_array(i);
i=i+1; p.NN         = data_array(i);
i=i+1; p.sigma      = data_array(i);
i=i+1; p.beta       = data_array(i);
i=i+1; p.alpha      = data_array(i);
i=i+1; p.delta      = data_array(i);
i=i+1; p.tauk       = data_array(i);
i=i+1; p.taul       = data_array(i);
i=i+1; p.theta      = data_array(i);
i=i+1; p.dummy      = data_array(i);    % skip rho from file
i=i+1; p.gamma      = data_array(i);
i=i+1; p.var_e      = data_array(i);
i=i+1; p.var_y1     = data_array(i);
i=i+1; p.n_pop      = data_array(i);
i=i+1; p.TFP        = data_array(i);
i=i+1; p.z_grid_opt = data_array(i);
i=i+1; p.z_min_w    = data_array(i);
i=i+1; p.z_max_w1   = data_array(i);
i=i+1; p.z_max_w2   = data_array(i);
i=i+1; p.a_min      = data_array(i);
i=i+1; p.a_max      = data_array(i);
i=i+1; p.a_scale    = data_array(i);
i=i+1; p.K0         = data_array(i);
i=i+1; p.T0         = data_array(i);
i=i+1; p.run_opt    = data_array(i);

% fclose returns 0 if successful, -1 if not
fclose(fid);

% derived parameters
p.sd_e = sqrt(p.var_e);
p.sd_y1 = sqrt(p.var_y1);

k = 4;
p.n_a2 = p.n_a + (p.n_a - 1)*k;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ y_bar_vals, s_next_vals ] = load_other_inputs( filename, n_age )

% fopen returns fid >= 0 if successful, -1 if not
fid = fopen(filename, 'r');

% read data, assume file is valid
data_cell = textscan(fid, '%d %f %f');

% fclose returns 0 if successful, -1 if not
fclose(fid);

age_vals    = data_cell{1};   % should be 1:n_age
y_bar_vals  = data_cell{2};
s_next_vals = data_cell{3};

% error checking
n1 = length(data_cell{1});
n2 = length(data_cell{2});
n3 = length(data_cell{3});

if (n1 ~= n_age || n2 ~= n_age || n3 ~= n_age)
    error('number of elements not equal to n_age')
end

if (s_next_vals(n_age) ~= 0)
    error('s(I+1) should equal zero')
end

if any( s_next_vals < 0 )
    error('s(i+1) should be greater than or equal to zero')
elseif( any(s_next_vals > 1) )
    error('s(i+1) should be less than or equal to 1')
end

end
