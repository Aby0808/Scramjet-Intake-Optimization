% This is the main script for the HCM intake optimization program
% The optimization alogirthm is Non-Dominated Sorting Genetic(NSGA)-II

% Yusoff, Yusliza, Mohd Salihin Ngadiman, and Azlan Mohd Zain. "Overview 
% of NSGA-II for optimizing machining process parameters." Procedia Engineering 15 (2011): 3978-3983.
% https://www.mathworks.com/help/gads/gamultiobj-algorithm.html

% This algorithm optimizes the first to fourth ramp deflection angles for min Drag and max Pressure Recovery
% Hard constraints are put on max length, min throat temperature and required mass flow rate

clc
clear
close all

projectPaths = setup_project_paths();

clear get_BL_SHAPE_PARAM get_FRSTM_TH_PARAM  % clear all getters
clear functions                              % clear all function caches

%% HCM intake: NSGA-II optimization (parallel, globals kept)

% ---------- pool setup ----------
delete(gcp('nocreate'));                     % clean any old pool
pool = parpool(8);                           % start fresh pool

% Make sure workers can see the helper and data
addAttachedFiles(pool, { ...
    fullfile(projectPaths.optimizationDir, 'FUN_setup_globals_fast.m'), ...  % helper to set globals on workers
    fullfile(projectPaths.optimizationDir, 'FUN_set_freestream_throat_params.m'), ...
    fullfile(projectPaths.thermoDir, 'FUN_get_prop_NASA9.m'), ...
    fullfile(projectPaths.optimizationDir, 'FUN_evaluateObjectives.m'), ...
    fullfile(projectPaths.dataDir, 'bl_shap_param.dat') });

% ---------- build params ON CLIENT (fast, cached) ----------
BL = get_BL_SHAPE_PARAM();                   % reads .dat once per session
FR = get_FRSTM_TH_PARAM();                   % computed once; cached

FUN_setup_globals_fast(BL, FR);

% ---------- initialize globals on all workers ----------
% (Workers cannot call local functions from a script, so we call the helper)
F = parfevalOnAll(@FUN_setup_globals_fast, 0, BL, FR);
wait(F);
for k = 1:numel(F)                          % surface any worker error
    if ~isempty(F(k).Error), rethrow(F(k).Error); end
end

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
testx = repmat(lb,5,1) + rand(5,4).*repmat(ub-lb,5,1);
futs  = parallel.FevalFuture.empty;
for i = 1:size(testx,1)
    futs(i) = parfeval(@FUN_evaluateObjectives, 1, testx(i,:));  %#ok<AGROW>
end
for i = 1:numel(futs)
    try
        fetchOutputs(futs(i));
    catch ME
        error("Worker preflight failed:\n%s", getReport(ME,'extended'));
    end
end

% ---------- run NSGA-II ----------
[x, fval] = gamultiobj(@FUN_evaluateObjectives, 4, [], [], [], [], lb, ub, options);

% ---------- persist results + params for post-processing ----------
save(fullfile(projectPaths.resultsDir, 'optimization_results.mat'), 'x', 'fval');

params.BL = BL;
params.FR = FR;
save(fullfile(projectPaths.resultsDir, 'run_params.mat'),'-struct','params');   % tiny MAT to reuse same freestream/throat

% ---------- cleanup ----------
delete(gcp('nocreate'));
disp('Optimization completed successfully. Parallel processing safely closed.');

%% ====== LOCAL GETTERS (cached; read/compute once per session) ======
function BL = get_BL_SHAPE_PARAM()
% Cached loader for boundary-layer shape parameters (worker-safe on client)
    persistent BL_cached
    if isempty(BL_cached)
        paths = setup_project_paths();
        dat  = fullfile(paths.dataDir,'bl_shap_param.dat');
        BL_cached = table2array(readtable(dat));
    end
    BL = BL_cached;
end

function FR = get_FRSTM_TH_PARAM()
% Cached builder for freestream/throat parameter vector/struct
    persistent FR_cached
    if isempty(FR_cached)
        FR_cached = FUN_set_freestream_throat_params();
    end
    FR = FR_cached;
end
