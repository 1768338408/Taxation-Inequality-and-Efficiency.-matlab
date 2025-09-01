% CALC_DIST_STATS  Compute macro and distributional statistics.
%
%   This function has been extended to account for two distinct lump‑sum transfers:
%   the accidental bequest transfer T_b and the government transfer Ts.  In
%   disposable income and other measures, the total transfer received by
%   households is T_b + Ts.
%
function [ stats, Phi_all ] = calc_dist_stats(a_grid, a_grid2, z_grid, y_bar_vals, Phi, mu,  ...
    w, r, K, T_b, Ts, L, ss_ben, p)

% Expect 14 input arguments: see above for descriptions.
if (nargin ~= 14)
    error('calc_stats: usage')
end

n_a2  = p.n_a2;
n_a   = p.n_a;
n_z   = p.n_z;
n_age = p.n_age;

stats = {};

% other variables
Y = p.TFP * (K^p.alpha) * (L^(1-p.alpha));  % output
G = p.tauk*(r*K) + p.taul*(w*L);                      % government spending

stats.K      = K;       % aggregate capital
% Save transfers separately: T_b is the accidental bequest transfer, Ts the government transfer
stats.T_b    = T_b;     % per‑capita accidental bequest transfer
stats.Ts     = Ts;      % per‑capita government transfer
stats.L      = L;       % aggregate labor
stats.Y      = Y;       % aggregate output
stats.G      = G;       % government spending
stats.K_Y    = K/Y;     % capital to output ratio
stats.K_L    = K/L;     % capital to labor ratio
stats.ss_ben = ss_ben;  % social security benefit
stats.w      = w;       % wage
stats.r      = r;       % interest rate

% compute distribution over assets and shock, collapsing across age
Phi_all = zeros(n_a2,n_z);
for i_a = 1:n_a2
    for i_z = 1:n_z
        Phi_all(i_a,i_z) = sum( squeeze(Phi(i_a,i_z,:)) .* mu );
    end
end

% approximate distribution on a grid (asset_grid) and compute the 
% PDF of households on this grid.
asset_grid = [a_grid2(1); 0.5*a_grid2(2:end)+0.5*a_grid2(1:end-1)];
dist_pdf = [Phi_all(1,:); Phi_all(2:end,:) - Phi_all(1:end-1,:)];

% construct variables needed to compute lorenz curve
dist_pdf_all = sum(dist_pdf, 2);

% compute Lorenz curve for wealth
[x, L] = calc_lorenz_curve(asset_grid, dist_pdf_all);

% wealth gini
stats.wealth_gini = calc_gini(x, L);

% fraction of wealth in top 1%, 5%, 20%, 40%, 60%, 80%
stats.wealth_top_01 = 1 - interp1(x, L, 1 - 0.01);
stats.wealth_top_05 = 1 - interp1(x, L, 1 - 0.05);
stats.wealth_top_20 = 1 - interp1(x, L, 1 - 0.20);
stats.wealth_top_40 = 1 - interp1(x, L, 1 - 0.40);
stats.wealth_top_60 = 1 - interp1(x, L, 1 - 0.60);
stats.wealth_top_80 = 1 - interp1(x, L, 1 - 0.80);

% fraction of wealth in bottom 20%, 40%
stats.wealth_bottom_20 = interp1(x, L, 0.20);
stats.wealth_bottom_40 = interp1(x, L, 0.40);

% share of households with zero wealth
stats.wealth_zero = interp1(a_grid2, sum(Phi_all,2), 0);

% earnings policy function
earn_pol = zeros(n_a, n_z, n_age);
for age = 1:n_age
    e_vals = y_bar_vals(age) * exp( z_grid );
    earn_pol(:,:,age) = repmat(w*e_vals', n_a, 1);
end

% compute earnings statistics
[dist_stats] = calc_policy_dist_stats(a_grid, a_grid2, Phi, mu, earn_pol, p);

% extract relevant statistics
stats.earnings_gini = dist_stats.pol_gini;
stats.earnings_top_01 = dist_stats.pol_top_01;
stats.earnings_top_05 = dist_stats.pol_top_05;
stats.earnings_top_20 = dist_stats.pol_top_20;
stats.earnings_top_40 = dist_stats.pol_top_40;
stats.earnings_top_60 = dist_stats.pol_top_60;
stats.earnings_top_80 = dist_stats.pol_top_80;
stats.earnings_bottom_20 = dist_stats.pol_bottom_20;
stats.earnings_bottom_40 = dist_stats.pol_bottom_40;
stats.earnings_gini_by_age = dist_stats.pol_gini_by_age;

% =====================================================================
% Compute wealth Gini *by age*
%
% To capture the distribution of assets within each age cohort, we treat
% the asset holdings as a policy variable and compute the Gini
% coefficients within each age group.  We construct a policy function
% where, for every productivity state and age, the policy equals the
% asset grid itself.  The helper function calc_policy_dist_stats
% computes the Lorenz curve and Gini coefficients by age.

% Construct a policy function for assets on the coarse asset grid.  The
% dimensions are n_a-by-n_z-by-n_age.  For each age and productivity
% state, the policy simply equals the asset grid.  This yields a
% distribution over assets identical to the underlying asset holdings.
asset_pol = zeros(p.n_a, p.n_z, p.n_age);
for age_idx = 1:p.n_age
    for iz_idx = 1:p.n_z
        asset_pol(:,iz_idx,age_idx) = a_grid(:);
    end
end

% Use the same routine to compute wealth distribution statistics.  This
% returns an object with overall and age-specific Gini coefficients.
asset_stats = calc_policy_dist_stats(a_grid, a_grid2, Phi, mu, asset_pol, p);

% Save the Gini-by-age vector for wealth.
stats.wealth_gini_by_age = asset_stats.pol_gini_by_age;

% =====================================================================
% Compute disposable income statistics
%
% Disposable income is net-of-tax interest income plus net-of-tax labor
% income plus lump-sum transfers and pension benefits.  For each asset
% level a and productivity state z, it is defined as:
%   y_disp = (1 - p.tauk) * r * a + (1 - p.theta - p.taul) * w * e(z,age)
%            + T + b(age)
% where b(age) is the social-security benefit for retirees (zero for
% working ages).  We build a policy matrix disp_pol(a,z,age) on the
% coarse grid, interpolate it to the fine grid, and compute the Lorenz
% curve via calc_policy_dist_stats, analogously to the earnings measure.

% initialise disposable income policy function
disp_pol = zeros(n_a, n_z, n_age);

% asset-related component (net interest) does not depend on z or age
inc_asset_part = (1 - p.tauk) * r * a_grid(:);  % n_a-by-1

for age = 1:n_age
    % labour efficiency for this age across z states
    % labour efficiency for this age across z states.  Ensure e_vals is a row
    % vector regardless of z_grid orientation by reshaping the exponential term.
    e_vals = (y_bar_vals(age) * exp( z_grid(:) ) )';   % 1-by-n_z
    % labour-related component (net wage). inc_labour_part will be 1-by-n_z
    inc_labour_part = (1 - p.theta - p.taul) * w * e_vals;
    % social security benefit: non‑zero only for retired
    if age < p.R_age
        b_age = 0.0;
    else
        b_age = ss_ben;
    end
    % total disposable income for each a,z
    % Note: use repmat to align dimensions: inc_asset_part (n_a-by-1) +
    % inc_labour_part (1-by-n_z) + scalar (T + b_age)
    disp_pol(:,:,age) = repmat(inc_asset_part, 1, n_z) + ...
                        repmat(inc_labour_part, n_a, 1) + ...
                        (T_b + Ts + b_age);
end

% compute disposable income distribution statistics
disp_stats = calc_policy_dist_stats(a_grid, a_grid2, Phi, mu, disp_pol, p);

% save relevant disposable income statistics
stats.disp_income_gini = disp_stats.pol_gini;
stats.disp_income_top_01 = disp_stats.pol_top_01;
stats.disp_income_top_05 = disp_stats.pol_top_05;
stats.disp_income_top_20 = disp_stats.pol_top_20;
stats.disp_income_top_40 = disp_stats.pol_top_40;
stats.disp_income_top_60 = disp_stats.pol_top_60;
stats.disp_income_top_80 = disp_stats.pol_top_80;
stats.disp_income_bottom_20 = disp_stats.pol_bottom_20;
stats.disp_income_bottom_40 = disp_stats.pol_bottom_40;
stats.disp_income_gini_by_age = disp_stats.pol_gini_by_age;

end

% -------------------------------------------------------------------------
% calc_policy_dist_stats
% -------------------------------------------------------------------------
function [dist_stats] = calc_policy_dist_stats(a_grid, a_grid2, Phi, mu, pol, p)

n_a2  = p.n_a2;
n_z   = p.n_z;
n_age = p.n_age;

pol2 = zeros(n_a2,n_z,n_age);

dist_stats = {};

% first interpolate policy function on finer grid
for age = 1:n_age
    for i_z = 1:n_z
        pol2(:,i_z,age) = interp1(a_grid, pol(:,i_z,age), a_grid2, 'linear');
    end
end

% approximate pdf/mass of agents at each grid point
dist_pdf = zeros(n_a2,n_z,n_age);
for age = 1:n_age
    dist_pdf(1,:,age) = Phi(1,:,age) * mu(age);
    dist_pdf(2:end,:,age) = (Phi(2:end,:,age) - Phi(1:end-1,:,age)) .* mu(age);
end

% construct overall Lorenz curve for this policy variable
[x_all, L_all] = calc_lorenz_curve_from_pdf_and_policy(pol2, dist_pdf);

% compute statistics from Lorenz curve
dist_stats.pol_gini = calc_gini(x_all, L_all);

dist_stats.pol_top_01 = 1 - interp1(x_all, L_all, 1 - 0.01);
dist_stats.pol_top_05 = 1 - interp1(x_all, L_all, 1 - 0.05);
dist_stats.pol_top_20 = 1 - interp1(x_all, L_all, 1 - 0.20);
dist_stats.pol_top_40 = 1 - interp1(x_all, L_all, 1 - 0.40);
dist_stats.pol_top_60 = 1 - interp1(x_all, L_all, 1 - 0.60);
dist_stats.pol_top_80 = 1 - interp1(x_all, L_all, 1 - 0.80);
    
% fraction of wealth in bottom 20%, 40%
dist_stats.pol_bottom_20 = interp1(x_all, L_all, 0.20);
dist_stats.pol_bottom_40 = interp1(x_all, L_all, 0.40);

% Now calculate gini within each age group
gini_vec = zeros(n_age,1);
for age = 1:n_age
    % pdf just for this age
    dist_pdf_age = dist_pdf(:,:,age);
    dist_pdf_age = dist_pdf_age ./ sum(dist_pdf_age(:));

    % policy just for this age
    pol2_age = pol2(:,:,age);

    % construct lorenz curve for this policy variable, within this age
    [x, L] = calc_lorenz_curve_from_pdf_and_policy(pol2_age, dist_pdf_age);

    % compute gini
    gini_vec(age) = calc_gini(x, L);

end

% save the gini for each age
dist_stats.pol_gini_by_age = gini_vec;


end

% -------------------------------------------------------------------------
% calc_lorenz_curve_from_pdf_and_policy
% -------------------------------------------------------------------------
function [x, L] = calc_lorenz_curve_from_pdf_and_policy(pol2, dist_pdf)
% pol2 is the policy function, defined on a finer grid for assets
% dist_pdf is the mass of agents at each grid point

% Re-shape policies and pdf distribution into a vector
pol_vec = pol2(:);
pdf_vec = dist_pdf(:);

% sort by policy value
[pol_sorted, idx] = sort(pol_vec);
pdf_sorted = pdf_vec(idx);

% compute Lorenz curve
[x, L] = calc_lorenz_curve(pol_sorted, pdf_sorted);

end

% -------------------------------------------------------------------------
% calc_lorenz_curve
% -------------------------------------------------------------------------
function [x, L] = calc_lorenz_curve(asset_grid, pdf)

asset_dist = asset_grid .* pdf;

A = sum(asset_dist);    % total assets (or other variables like income)
M = sum(pdf);           % total mass of agents

% cumulative share of assets
L = [0; cumsum(asset_dist) ./ A];
% cumulative share of agents
x = [0; cumsum(pdf) ./ M];

% remove duplicates
[x, idx_u] = unique(x);
L = L(idx_u);
end

% -------------------------------------------------------------------------
% calc_gini
% -------------------------------------------------------------------------
function [gini] = calc_gini(x, L)

dbg_flag = 0;
if (dbg_flag > 0)
    plot(x, L, 'b*-');
    grid on
end

pp_L = interp1(x, L, 'linear', 'pp');

L_fun = @(x) ppval(pp_L, x);

fun = @(x) x - L_fun(x);

%gini = 2*quad(fun, 0, 1);
gini = 2*integral(fun, 0, 1);

% another formula
%gini2 = 1 - 2*sum(L)/length(L);
    
end
