% CALC_DIST_HUGGETT  Simulate the distribution of agents and compute
% aggregate capital and accidental bequest transfers.
%
% The household budget includes both a government transfer (Ts) and an
% accidental bequest transfer (T_b).  Ts is passed in exogenously,
% T_b_val is the current iteration's accidental bequest transfer.  At
% the end of the simulation a new T_b_val_new is computed to update
% the fixed‑point.
%
function [K_val, T_b_val_new, avgs_by_age, Phi] = calc_dist_huggett(a_grid, z_grid, a_grid2, Pz_cdf, Pz1_cdf, ...
        mu, y_bar_vals, s_next_vals, g, w, r, ss_ben, Ts_val, T_b_val, p)
    
N     = p.NN; % number of consumers, per age group
n_z   = p.n_z;
n_age = p.n_age;
n_a2  = p.n_a2;

% set seed
rng(1985)

avgs_by_age = {};
avgs_by_age.cons         = zeros(n_age,1);  % consumption
avgs_by_age.earn         = zeros(n_age,1);  % earnings, w*e
avgs_by_age.K            = zeros(n_age,1);  % capital/savings, initial
avgs_by_age.Kn           = zeros(n_age,1);  % capital/savings, next-period
avgs_by_age.inc          = zeros(n_age,1);  % total pre-tax income
avgs_by_age.acc_bequests = zeros(n_age,1);  % accidental bequests

% simulated values for each age group
sim_an = zeros(N,1);
sim_y  = zeros(N,1);
sim_c  = zeros(N,1);

% vectors for additional simulated values
sim_e = zeros(N,1);

% CDF of distribution
Phi = zeros(n_a2, n_z, n_age);

% probability of surviving to current age
% s_this_vals(age) == probability survive to this age
s_this_vals = zeros(n_age,1);
s_this_vals(1) = 1;
s_this_vals(2:n_age) = s_next_vals(1:n_age-1);

% initialize values for simulated agents
sim_a = zeros(N,1); % age = 1 agents start with no assets
u_vec = rand([N,1]);
%sim_z = simulate_initial_markov_shocks(N, Pz1_cdf, u_vec);
sim_z = simulate_initial_markov_shocks_fast(N, n_z, Pz1_cdf, u_vec);

for age = 1:n_age
    % calculate average capital for this age group:
    avgs_by_age.K(age) = mean(sim_a);  % should be zero for age == 1

    % calculate accidental bequests
    avgs_by_age.acc_bequests(age) = (1-s_this_vals(age)) * (1+r*(1-p.tauk)) * mean(sim_a);

    % compute social security benefit for this age
    if (age < p.R_age)
        b_val = 0.0;
    else
        b_val = ss_ben;
    end

    % update distribution for next age group
    for i_z = 1:n_z
        % identify agents with productivity i_z
        idx = (sim_z == i_z);

        % determine asset choice for these agents
        sim_an(idx) = interp1(a_grid,  g(:,i_z,age), sim_a(idx));

        % compute labor endowment for this individual
        e_val = y_bar_vals(age) * exp(z_grid(i_z));
        sim_e(idx) = e_val;
        % compute total resources: net interest income + net wage income
        % plus government lump‑sum transfer Ts and accidental bequest transfer T_b
        sim_y(idx) = (1 + (1-p.tauk)*r)*sim_a(idx) + ...
                     (1-p.theta-p.taul)*w*e_val + ...
                     Ts_val + T_b_val + b_val;
        % compute consumption
        sim_c(idx) = sim_y(idx)  - sim_an(idx);
    end

    % additional averages for this age group
    avgs_by_age.cons(age) = mean(sim_c);
    avgs_by_age.earn(age) = mean(w*sim_e);
    avgs_by_age.inc(age)  = mean(r*sim_a + w*sim_e);
    avgs_by_age.Kn(age)   = mean(sim_an);

    % now simulate next period's z:        
    u_vec = rand([N,1]);
    % slower way:
    %sim_zn = simulate_markov_shocks(N, sim_z, Pz_cdf, u_vec);
    % faster way:
    sim_zn = simulate_markov_shocks_fast(n_z, sim_z, Pz_cdf, u_vec);

    % slow method to compute CDF of distribution
    %for i_z = 1:n_z        
    %    for i_a = 1:n_a2
    %        Phi(i_a,i_z,age) = sum(sim_a <= a_grid2(i_a) & sim_z == i_z)/N;
    %    end % i_a loop        
    %end % i_z loop
    % faster method to compute CDF of distribution
    for i_z = 1:n_z
        idx = (sim_z == i_z);
        sim_a_idx = sim_a(idx);
        if isempty(sim_a_idx)
            Phi(:,i_z,age) = 0;
        else
            % sim_a_idx' is a [1 x M] row vector, a_grid2 is a [n_a2 x 1] column vector
            % Matlab impliicly expands the expression inside the sum to a [n_a2 x M] matrix
            Phi(:,i_z,age) = sum(sim_a_idx' <= a_grid2, 2) ./ N;                
        end
    end

    if ( mod(age, 10) == 0 )
        fprintf('age = %d\n', age);
    end

    % update initial values
    sim_a = sim_an;
    sim_z = sim_zn;

end % age loop

% calculate total capital
K_val = sum( mu .* avgs_by_age.K );

% calculate per‑capita accidental bequests (transfer) to be used as T_b in next iteration
T_b_val_new = sum( mu .* avgs_by_age.acc_bequests );

end

% -------------------------------------------------------------------------
% function simulate_initial_markov_shocks
% -------------------------------------------------------------------------
function [sim_z] = simulate_initial_markov_shocks(N, Pz1_cdf, u_vec)

sim_z = zeros(N,1);

for i = 1:N 
    sim_z(i) = calc_markov_shock_from_uniform_rv(u_vec(i), Pz1_cdf);
end

end

% -------------------------------------------------------------------------
% function simulate_markov_shocks
% -------------------------------------------------------------------------
function [sim_zn] = simulate_markov_shocks(N, sim_z, Pz_cdf, u_vec)

sim_zn = zeros(N,1);
%u_vec = rand([N,1]);
for i = 1:N 
    sim_zn(i) = calc_markov_shock_from_uniform_rv(u_vec(i), Pz_cdf(sim_z(i),:));
end

end

% -------------------------------------------------------------------------
% function simulate_initial_markov_shocks
% -------------------------------------------------------------------------
function [sim_z] = simulate_initial_markov_shocks_fast(N, n_z, Pz1_cdf, u_vec)

U_vec = kron(u_vec, ones(1,n_z));
sim_z = sum( U_vec >= ones(N,1)*Pz1_cdf', 2 ) + 1;

end

% -------------------------------------------------------------------------
% function simulate_markov_shocks
% -------------------------------------------------------------------------
function [sim_zn] = simulate_markov_shocks_fast(n_z, sim_z, Pz_cdf, u_vec)

U_vec = kron(u_vec, ones(1,n_z));
sim_zn = sum( U_vec >= Pz_cdf(sim_z,:), 2 ) + 1;

end


% -------------------------------------------------------------------------
% function calc_markov_shock_from_uniform_rv
% -------------------------------------------------------------------------
function [ i_zn ] = calc_markov_shock_from_uniform_rv(u, Pz_cdf)

i_zn = 0;
done = false;
    
while (~done)
    
    % get cumulative probability
    i_zn = i_zn + 1;
   
    % stop if u < p
    % u has to be less than 1
    % at end, p = 1, so eventually this will fail
    done = (u < Pz_cdf(i_zn));
    
end
    
end

