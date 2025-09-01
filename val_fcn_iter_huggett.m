%
% Modified to include separate government transfer (Ts) and accidental bequest transfer (T_b).
% The household budget constraint now reads:
%   c + a' = a[1 + r(1-τ_k)] + (1 - θ - τ_ℓ)w e_j(z) + Ts + T_b + b_j
% where Ts is the per‑capita government transfer financed by taxes and T_b
% is the per‑capita transfer of accidental bequests.  Both are passed
% explicitly as arguments.
%
function [ v, g ] = val_fcn_iter_huggett( a_grid, z_grid, Pz, ...
        y_bar_vals, s_next_vals, Ts, T_b, w, r, p )

% ------------------------------------------------------------------------
% Value-function iteration for Huggett life-cycle model
% (Retirement benefit = ρ · last-period wage, state-specific)
% ------------------------------------------------------------------------
n_a   = p.n_a;
n_z   = p.n_z;
n_age = p.n_age;

v  = zeros(n_a,n_z,n_age);
g  = zeros(n_a,n_z,n_age);
v0 = zeros(n_a,n_z);          % next-period value

dbg_flag = 0;

tauk  = p.tauk;     % capital-income tax
taul  = p.taul;     % labor-income  tax
theta = p.theta;    % social-security contribution rate

for age = n_age:-1:1
    
    % probability of surviving to age+1
    s_next = s_next_vals(age);

    % continuation value  V(a',z')  (same for all i_z)
    cv = v0 * Pz';
    
    % -------- loop over productivity states ----------------------------
    for i_z = 1:n_z
        
        % retirement benefit: 0 if working, otherwise ρ·last wage
        if age < p.R_age
            b_val = 0.0;
        else
            e_last = y_bar_vals(p.R_age-1) * exp( z_grid(i_z) );
            b_val  = p.rho * w * e_last;   
        end
        
        % current labour endowment
        e_val = y_bar_vals(age) * exp( z_grid(i_z) );
        
        % total resources (rows index current a)
        % Household receives net interest income, net wage income,
        % the government lump‑sum transfer Ts, the accidental bequest transfer T_b,
        % and any retirement benefit b_val.
        y_grid = (1 + (1-tauk)*r)*a_grid ...
               + (1-theta-taul)*w*e_val ...
               + Ts + T_b + b_val;
        
        % consumption matrix  c(a',a)
        c_choices = ones(n_a,1)*y_grid' - a_grid*ones(1,n_a);
        
        % utility matrix
        if (p.sigma == 1)
            u_choices = log(c_choices);
        else
            u_choices = c_choices.^(1-p.sigma) / (1-p.sigma);
        end
        
        % value matrix
        V_choices = u_choices + p.beta * s_next * cv(:,i_z);
        V_choices(c_choices <= 0) = -Inf;     % infeasible
        
        % optimal choice for each a
        [vt, gt]     = max(V_choices);
        v(:,i_z,age) = vt';
        g(:,i_z,age) = a_grid(gt');
    end
    % -------------------------------------------------------------------
    
    % optional debug output
    if dbg_flag
        if mod(age,10)==0, fprintf('age = %d\n',age); end
    else
        if mod(age,10)==0, fprintf('age = %d\n',age); end
    end
    
    v0 = v(:,:,age);    % update V_{t+1}
end
end

