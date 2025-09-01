function mu = calc_age_distribution(p, s_next_vals)
% CALC_AGE_DISTRIBUTION  Compute stationary age weights.
%
%   mu = CALC_AGE_DISTRIBUTION(p, s_next_vals) returns a column vector mu
%   whose i-th element is the fraction of the population that is exactly
%   age i, given the survival probabilities in s_next_vals.
%
%   Inputs
%   ------
%   p             : struct, must contain the field n_age (scalar)
%   s_next_vals   : n_age-by-1 vector where s_next_vals(i) equals the
%                   conditional probability of surviving from age i to
%                   age i+1.  By construction s_next_vals(n_age) = 0.
%
%   Output
%   ------
%   mu            : n_age-by-1 vector with sum(mu) = 1.
%

% Number of discrete ages
n_age = p.n_age;

% 1.  Un-normalised mass: set age-1 mass to 1, then propagate forward
mu = ones(n_age, 1);
for age = 2:n_age
    mu(age) = mu(age-1) * s_next_vals(age-1);
end

% 2.  Normalise so that the masses sum to one
mu = mu / sum(mu);

end
