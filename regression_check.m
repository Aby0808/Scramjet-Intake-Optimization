clear
clc

projectPaths = setup_project_paths();
baselineFile = fullfile(projectPaths.resultsDir, 'regression_baseline.mat');

BL = table2array(readtable(fullfile(projectPaths.dataDir, "bl_shap_param.dat")));
FR = FUN_set_freestream_throat_params();
FUN_setup_globals_fast(BL, FR);

cases = buildRegressionCases();
results = evaluateCases(cases, FR);

if isfile(baselineFile)
    baselineData = load(baselineFile, 'results', 'metadata');
    compareResults(baselineData.results, results, baselineData.metadata);
else
    metadata.createdOn = char(datetime('now'));
    metadata.caseLabels = string({results.label});
    metadata.notes = "Baseline generated from the current trusted implementation.";
    save(baselineFile, 'results', 'metadata');
    fprintf('Saved baseline to %s\n', baselineFile);
    fprintf('Run this script again after refactoring to compare against the baseline.\n');
end

function cases = buildRegressionCases()
cases = struct('label', {}, 'theta', {});

optimizedCases = loadOptimizedCases();
for i = 1:numel(optimizedCases)
    cases(end+1).label = optimizedCases(i).label; %#ok<SAGROW>
    cases(end).theta = optimizedCases(i).theta;
end

manualThetas = [
    0.1, 0.1, 0.1, 0.1
    6.0, 4.0, 3.0, 1.0
    25.0, 25.0, 25.0, 25.0
];

manualLabels = ["manual_low_angles", "manual_mid_angles", "manual_high_angles"];

for i = 1:size(manualThetas, 1)
    cases(end+1).label = manualLabels(i); %#ok<SAGROW>
    cases(end).theta = manualThetas(i, :);
end
end

function optimizedCases = loadOptimizedCases()
optimizedCases = struct('label', {}, 'theta', {});
projectPaths = setup_project_paths();
candidateFiles = [ ...
    fullfile(projectPaths.resultsDir, "optimization_results_w_int_constraints.mat"), ...
    fullfile(projectPaths.resultsDir, "optimization_results.mat"), ...
    fullfile(projectPaths.resultsDir, "optimization_results_wo_int_constraints.mat"), ...
    fullfile(projectPaths.resultsDir, "optimization_results_no_int_constrain.mat")];

for i = 1:numel(candidateFiles)
    fileName = candidateFiles(i);
    if ~isfile(fileName)
        continue
    end

    data = load(fileName);
    if ~isfield(data, 'x') || isempty(data.x)
        continue
    end

    indices = unique([1, max(1, round(size(data.x,1)/2)), size(data.x,1)]);
    labels = ["opt_first", "opt_middle", "opt_last"];
    for j = 1:numel(indices)
        optimizedCases(end+1).label = labels(j); %#ok<SAGROW>
        optimizedCases(end).theta = data.x(indices(j), :);
    end
    return
end
end

function results = evaluateCases(cases, FR)
alpha = FR.alpha;
M_oo = FR.M_oo;
P_oo = FR.P_oo;
T_oo = FR.T_oo;
m_dot = FR.m_dot;
h_th = FR.h_th;
T_th = FR.T_th;

penalty = [1e5, -1e5];
results = repmat(struct( ...
    'label', "", ...
    'theta', zeros(1,4), ...
    'objectiveFunction', zeros(1,2), ...
    'optimizerObjective', zeros(1,2), ...
    'postprocess', zeros(1,6), ...
    'isPenalty', false), 1, numel(cases));

fprintf('Evaluating %d regression cases...\n', numel(cases));

for i = 1:numel(cases)
    theta = cases(i).theta;
    [objFn, ~, ~, postprocess] = FUN_objective_function(alpha, M_oo, P_oo, T_oo, ...
        m_dot, h_th, T_th, theta, 'n', 15);
    optimizerObjective = FUN_evaluateObjectives(theta);

    results(i).label = string(cases(i).label);
    results(i).theta = theta;
    results(i).objectiveFunction = objFn;
    results(i).optimizerObjective = optimizerObjective;
    results(i).postprocess = postprocess;
    results(i).isPenalty = isequal(size(objFn), size(penalty)) && all(abs(objFn - penalty) < 1e-12);

    fprintf('  %-18s theta = [%6.3f %6.3f %6.3f %6.3f]\n', results(i).label, theta);
end
end

function compareResults(baseline, current, metadata)
if numel(baseline) ~= numel(current)
    error('Regression case count changed: baseline has %d cases, current run has %d.', ...
        numel(baseline), numel(current));
end

objTol = 1e-8;
metricTol = 1e-8;
allPass = true;

fprintf('Comparing against baseline created on %s\n', metadata.createdOn);

for i = 1:numel(current)
    if baseline(i).label ~= current(i).label
        error('Case order changed at index %d: baseline=%s, current=%s.', ...
            i, baseline(i).label, current(i).label);
    end

    if any(abs(baseline(i).theta - current(i).theta) > 1e-12)
        error('Theta changed for case %s.', current(i).label);
    end

    if baseline(i).isPenalty ~= current(i).isPenalty
        allPass = false;
        fprintf('FAIL %-18s penalty flag changed.\n', current(i).label);
        continue
    end

    objDiff = max(abs(baseline(i).objectiveFunction - current(i).objectiveFunction));
    evalDiff = max(abs(baseline(i).optimizerObjective - current(i).optimizerObjective));
    postDiff = max(abs(baseline(i).postprocess - current(i).postprocess));

    if objDiff > objTol || evalDiff > objTol || postDiff > metricTol
        allPass = false;
        fprintf(['FAIL %-18s objDiff=% .3e evalDiff=% .3e postDiff=% .3e\n'], ...
            current(i).label, objDiff, evalDiff, postDiff);
    else
        fprintf('PASS %-18s\n', current(i).label);
    end
end

if allPass
    fprintf('All regression cases matched the saved baseline.\n');
else
    error('Regression mismatch detected. Review the failing cases above.');
end
end
