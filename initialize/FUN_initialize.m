function initData = FUN_initialize()
% Bundle the driver-side initialization required before calling gamultiobj.
% This collects path setup, cached data loading, pool startup, file
% attachment, and global initialization into one place so the main script
% reads like the optimization workflow rather than setup plumbing.

projectPaths = setup_project_paths();

clear FUN_get_BL_SHAPE_PARAM FUN_get_FRSTM_TH_PARAM  % clear all getters
clear functions                                      % clear all function caches

pool = FUN_start_parallel_pool(8);

% Make sure workers can see the helper and data
addAttachedFiles(pool, { ...
    fullfile(projectPaths.initializeDir, 'FUN_setup_globals_fast.m'), ...  % helper to set globals on workers
    fullfile(projectPaths.initializeDir, 'FUN_set_freestream_throat_params.m'), ...
    fullfile(projectPaths.thermoDir, 'FUN_get_prop_NASA9.m'), ...
    fullfile(projectPaths.optimizationDir, 'FUN_evaluateObjectives.m'), ...
    fullfile(projectPaths.dataDir, 'bl_shap_param.dat') });

% Build params ON CLIENT (fast, cached)
BL = FUN_get_BL_SHAPE_PARAM();   % reads .dat once per session
FR = FUN_get_FRSTM_TH_PARAM();   % computed once; cached

FUN_setup_globals_fast(BL, FR);

% Initialize globals on all workers.
% Workers cannot call local functions from a script, so we call the helper.
F = parfevalOnAll(@FUN_setup_globals_fast, 0, BL, FR);
wait(F);
for k = 1:numel(F)  % surface any worker error
    if ~isempty(F(k).Error), rethrow(F(k).Error); end
end

initData = struct( ...
    'projectPaths', projectPaths, ...
    'pool', pool, ...
    'BL', BL, ...
    'FR', FR);
end

