%% run_calibration.m  –  one-click driver (SMLT version)
% ------------------------------------------------------------
% 0.  Clean workspace
% 1.  Check Statistics & ML Toolbox (normcdf, lsqnonlin)
% 2.  Delegate to main.m  (main controls solve/calibrate/plot)
% ------------------------------------------------------------
clearvars; close all; clc;

%% 1) Quick sanity check: key SMLT functions
funcs = {'normcdf','lsqnonlin','optimoptions'};
for k = 1:numel(funcs)
    assert(exist(funcs{k},'file') == 2, ...
        'Function %s not found. Check Statistics and Machine Learning Toolbox.', funcs{k});
end
disp('✓ Statistics & ML Toolbox detected – proceeding');

%% 2) Launch the top-level script
if exist('main','file') ~= 2
    error('main.m not found in the current directory.');
end

fprintf('\n=========  Huggett model run begins  =========\n');
tic;
main;                          % main.m decides run_opt (1=SS, 2=calibration)
toc;
fprintf('=========  Huggett model run finished =========\n');
