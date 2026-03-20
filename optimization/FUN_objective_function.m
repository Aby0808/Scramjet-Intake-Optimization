function [obj_fn, ramp_coord, coord_ramp1, postprocess] = FUN_objective_function(alpha, M_oo, P_oo, ...
    T_oo, m_dot, h_th, T_t, theta, post, config_param)

% Compatibility wrapper around the full design-evaluation pipeline.
result = FUN_evaluate_design(alpha, M_oo, P_oo, T_oo, m_dot, h_th, T_t, theta, post, config_param);

obj_fn = result.objectiveFunction;
ramp_coord = result.rampCoord;
coord_ramp1 = result.coordRamp1;
postprocess = result.postprocess;

end
