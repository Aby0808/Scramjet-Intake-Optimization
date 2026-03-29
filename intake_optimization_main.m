% This is the main script for the HCM intake optimization program
% The optimization alogirthm is Non-Dominated Sorting Genetic(NSGA)-II

% This algorithm optimizes the first to fourth ramp deflection angles for min Drag and max Pressure Recovery
% Hard constraints are put on max length, min throat temperature and required mass flow rate

clc
clear
close all

setup_project_paths();
initData = FUN_initialize();

%% HCM intake: NSGA-II optimization (parallel, globals kept)

% ---------- bounds ----------
% set range for input parameters for optimization(ramp 1 to 4 angles)
% [theta1, theta2, theta3, theta4]  in degrees
lb = [0.1, 0.1, 0.1, 0.1];
ub = [25.0, 25.0, 25.0, 25.0];

% ---------- GA options ----------
rng(42);  % reproducible
options = optimoptions('gamultiobj', ...
    'PopulationSize',      300, ...
    'MaxGenerations',      300, ...
    'CrossoverFraction',   0.90, ...
    'FunctionTolerance',   1e-4, ...
    'ConstraintTolerance', 1e-6, ...
    'ParetoFraction',      0.60, ...
    'UseParallel',         true, ...
    'UseVectorized',       false, ...
    'Display',             'iter', ...
    'PlotFcn',             {@gaplotpareto, @gaplotparetodistance});

% ---------- smoke test (catches worker issues early) ----------
FUN_run_worker_smoke_test(lb, ub);

% ---------- run NSGA-II ----------
% Yusoff, Yusliza, Mohd Salihin Ngadiman, and Azlan Mohd Zain. "Overview 
% of NSGA-II for optimizing machining process parameters." Procedia Engineering 15 (2011): 3978-3983.
% https://www.mathworks.com/help/gads/gamultiobj-algorithm.html

[x, fval] = gamultiobj(@FUN_evaluateObjectives, 4, [], [], [], [], lb, ub, options);

% ---------- persist results + params for post-processing ----------
FUN_finalize(initData, x, fval);
