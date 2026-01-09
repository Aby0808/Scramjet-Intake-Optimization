clc
clear
close all

% Clear existing parallel pool (if any) to prevent conflicts
delete(gcp('nocreate')); 

% Start a parallel pool with 6 cores
parpool(7); 

% Set bounds for ramp angles
lb = [0.01, 0.01, 0.01, 0.01, 0.5];   % Minimum values for ramp angles
ub = [40, 40, 40, 40, 1]; % Maximum values for ramp angles

% Set NSGA-II options with parallel execution enabled
options = optimoptions('gamultiobj', ...
    'PopulationSize', 400, ...  % Small population for quick visualization
    'MaxGenerations', 200, ...  % Lower generations to start with
    'CrossoverFraction', 0.8, ...
    'Display', 'iter', ...
    'UseParallel', true); % Enable parallel processing

% Run NSGA-II optimization
[x, fval] = gamultiobj(@FUN_evaluateObjectives, 5, [], [], [], [], lb, ub, options);

% Save results
save('optimization_results.mat', 'x', 'fval');

% Close the parallel pool after computation ends (cleaning up resources)
delete(gcp('nocreate')); 

disp('Optimization completed successfully. Parallel processing safely closed.');
