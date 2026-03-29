function pool = FUN_start_parallel_pool(numWorkers)
% Start from a clean pool state so each optimization run is reproducible
% and worker paths/globals are initialized in a controlled way.

if nargin < 1
    numWorkers = 8;
end

delete(gcp('nocreate'));   % clean any old pool
pool = parpool(numWorkers); % start fresh pool
end
