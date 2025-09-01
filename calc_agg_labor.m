function L = calc_agg_labor(n_z, n_age, z_grid, y_bar_vals, mu, Pz1, Pz)
% CALC_AGG_LABOR   Aggregate effective labor supply.
%
%   L = CALC_AGG_LABOR(n_z, n_age, z_grid, y_bar_vals, mu, Pz1, Pz)
%
%   Computes total effective labor L by integrating over age and
%   productivity shocks.
%
%   Inputs
%     n_z         - number of productivity states
%     n_age       - number of discrete age groups
%     z_grid      - n_z-by-1 vector, log productivity grid
%     y_bar_vals  - n_age-by-1 vector, deterministic efficiency profile
%     mu          - n_age-by-1 population mass at each age (sum = 1)
%     Pz1         - n_z-by-1 initial distribution over z at age 1
%     Pz          - n_z-by-n_z Markov transition matrix of z
%
%   Output
%     L           - scalar, aggregate effective labor
%

% ---------- 1) Distribute productivity across ages ----------------------
dist_z          = zeros(n_z, n_age);
dist_z(:,1)     = Pz1(:);          % age-1 distribution
for a = 2:n_age
    dist_z(:,a) = Pz' * dist_z(:,a-1);
end

% ---------- 2) Expected exp(z) conditional on age -----------------------
Eexpz_by_age = dist_z' * exp(z_grid);   % n_age-by-1

% ---------- 3) Aggregate effective labor --------------------------------
L = sum( mu .* y_bar_vals .* Eexpz_by_age );

end
