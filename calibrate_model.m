function [x] = calibrate_model(p, y_bar_vals, s_next_vals)

% Set to 1 to load previously computed SMM solution
% Set to 0 otherwise
load_smm_solution = 0;

% data moments
moments_data = [ ...
    2.90;      % K/Y
    0.49;   % earnings gini    
    ];

% set weight matrix for SMM criterion
WgtMtx = diag(1./moments_data.^2);

% upper cholesky decomposition so that we can
% re-write SMM criterion as a sum of squares
WgtMtxU = chol(WgtMtx);

% initial guess for parameters and bounds:
A = [ ...
    % param   lb     ub
    p.beta,  0.92,  1.08;
    p.sigma,  1.0,   5.0; 
];


x0 = A(:,1);  % initial guess
lb = A(:,2);  % lower bound for parameters
ub = A(:,3);  % upper bound for parameters

% save parameter names
num_params = length(x0);
param_names = cell(num_params,1);
ip = 0;
ip = ip+1; param_names{ip} = 'beta';
ip = ip+1; param_names{ip} = 'sigma';

% save moment names
num_moments = length(moments_data);
moment_names = cell(num_moments,1);
im = 0;
im = im+1; moment_names{im} = 'K/Y';
im = im+1; moment_names{im} = 'earnings gini';

% Minimize SMM criterion to calibrate model
if (load_smm_solution > 0)
    load('smm_solution.mat', 'x', 'resnorm', 'moments_err', 'moments_model')
else
    fun = @(x) smm_criterion(x, p, y_bar_vals, s_next_vals, moments_data, WgtMtxU, param_names, moment_names);
    
    % Potentially choose different algorithm
    options = optimoptions(@lsqnonlin,...
        'Algorithm','trust-region-reflective', ... 
        'FunctionTolerance', 1e-4, ...
        'FiniteDifferenceStepSize', 0.01, ...
        'Display', 'iter' ...
        );
    %'Algorithm', 'levenberg-marquardt', ...

    % Minimize SMM criterion
    [x, resnorm] = lsqnonlin(fun, x0, lb, ub, options);

    % solve model again to determine get model moments
    fprintf('=============================\n')
    fprintf('SOLVING MODEL AGAIN:\n')
    fprintf('=============================\n')
    [moments_err, moments_model] = smm_criterion(x, p, y_bar_vals, s_next_vals, ...
        moments_data, WgtMtxU, param_names, moment_names);

    save('smm_solution.mat')
end

fprintf('=============================\n')
fprintf('SMM SOLUTION:\n')
fprintf('=============================\n')
print_parameters(x, param_names)
fprintf('\n')
print_moments(moments_model, moments_data, moment_names)
fprintf('\n')
fprintf('resnorm: %0.16e\n', resnorm)
%fprintf('resnorm: %0.16e\n', sum(moments_err.^2))

end

% -------------------------------------------------------------------------
% smm_criterion
% -------------------------------------------------------------------------
function [moments_err, moments_model] = smm_criterion(x, p, y_bar_vals, s_next_vals, ...
    moments_data, WgtMtxU, param_names, moment_names)

fprintf('+++++++++++++++++++++++++++++++++++++++++\n')
fprintf('Input Parameters:\n')
print_parameters(x, param_names)
fprintf('+++++++++++++++++++++++++++++++++++++++++\n')

% copy parameters
p_copy = p;

% set new parameters
p_copy.beta = x(1);
p_copy.sigma = x(2);

% solve model
results_opt = 0;
[stats] = solve_model(p_copy, y_bar_vals, s_next_vals, results_opt);

% get model moments
moments_model = [
    stats.K_Y;
    stats.earnings_gini
    ];

% compute error
moments_err = WgtMtxU * (moments_model - moments_data);

fprintf('+++++++++++++++++++++++++++++++++++++++++\n')
fprintf('Input Parameters:\n')
print_parameters(x, param_names)
fprintf('Solution:\n')
fprintf('K = %.16e\n', stats.K);
% Print both accidental and government lump‑sum transfers if available
if isfield(stats,'T_b')
    fprintf('T_b (accidental bequests) = %.16e\n', stats.T_b);
end
if isfield(stats,'Ts')
    fprintf('Ts (government transfer) = %.16e\n', stats.Ts);
end
fprintf('Model Moments:\n')
print_moments(moments_model, moments_data, moment_names)
fprintf('Residual Error: %.16e\n', sum(moments_err.^2))
fprintf('+++++++++++++++++++++++++++++++++++++++++\n')

end

% -------------------------------------------------------------------------
% print_parameters
% -------------------------------------------------------------------------
function [] = print_parameters(x, param_names)

num_params = length(x);

fprintf('%12s  %16s\n', 'PARAMETER', 'VALUE');
for ip = 1:num_params
    fprintf('%12s  %.16f\n', param_names{ip}, x(ip))
end

end

% -------------------------------------------------------------------------
% print_moments
% -------------------------------------------------------------------------
function [] = print_moments(moments_model, moments_data, moment_names)

num_moments = length(moments_model);

fprintf('%15s %10s %10s\n', 'MOMENT', 'MODEL', 'DATA')
for im = 1:num_moments
    fprintf('%15s %10.3f %10.3f\n', moment_names{im}, moments_model(im), moments_data(im));
end

end
