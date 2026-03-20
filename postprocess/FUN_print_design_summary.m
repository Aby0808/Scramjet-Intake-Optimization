function FUN_print_design_summary(P_oo, T_oo, P0_oo, ramp_param, Tth, Mth, D, L, lng_in, wdt_in, wdt_cwl, SI)
% Print the summary metrics for a post-processing run.
% This keeps reporting separate from the core design-evaluation logic.

fprintf('\nPR\t=\t%f', ramp_param(5,4)/P_oo)
% fprintf('\nTR\t=\t%f', ramp_param(5,5)/T_oo)
fprintf('\nTR\t=\t%f', Tth/T_oo)
fprintf('\nTPR\t=\t%f', ramp_param(5,6)/P0_oo)
fprintf('\nMth\t=\t%f', Mth)
fprintf('\nTth (K)\t=\t%f', Tth)
fprintf('\nDrag (N)\t=\t%f', D)
fprintf('\nL/D\t=\t%f', L/D)
fprintf('\nLength (m)\t=\t%f', lng_in)
fprintf('\nWidth intake (m)\t=\t%f', wdt_in)
fprintf('\nWidth cowl (m)\t=\t%f', wdt_cwl)
fprintf('\nStartabililty index\t=\t%f\n\n', SI)
end
