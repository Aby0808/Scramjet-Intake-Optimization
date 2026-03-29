function FUN_run_worker_smoke_test(lb, ub)
% Smoke test a few random points before launching the full optimization.
% This catches worker/path/global initialization issues early, when the
% failure is cheaper and easier to interpret.

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
end
