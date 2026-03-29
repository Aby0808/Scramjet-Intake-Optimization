function BL = FUN_get_BL_SHAPE_PARAM()
% Cached loader for boundary-layer shape parameters (worker-safe on client).
% Keeping this separate from the driver makes initialization reusable across
% optimization, regression checks, and any future batch workflows.

persistent BL_cached

if isempty(BL_cached)
    paths = setup_project_paths();
    dat  = fullfile(paths.dataDir, 'bl_shap_param.dat');
    BL_cached = table2array(readtable(dat));
end

BL = BL_cached;
end
