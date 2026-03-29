function FUN_finalize(initData, x, fval)
% Persist optimization outputs and close the parallel resources cleanly.
% Keeping this separate makes it easier to add more result packaging later
% without re-expanding the main driver script.

save(fullfile(initData.projectPaths.resultsDir, 'optimization_results.mat'), 'x', 'fval');

params.BL = initData.BL;
params.FR = initData.FR;
save(fullfile(initData.projectPaths.resultsDir, 'run_params.mat'), '-struct', 'params');   % tiny MAT to reuse same freestream/throat

delete(gcp('nocreate'));
disp('Optimization completed successfully. Parallel processing safely closed.');
end
