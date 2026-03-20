function obj = FUN_evaluateObjectives(param)
theta = param;  %(1:end-1);

%% ---- fixed parameters ------
% freestream and some fixed throat requirements

global FRSTM_TH_PARAM BL_SHAPE_PARAM

if isempty(BL_SHAPE_PARAM) || isempty(FRSTM_TH_PARAM)
    % incase globals are not set 
    BL = get_BL_SHAPE_PARAM();
    FR = get_FRSTM_TH_PARAM();
    FUN_setup_globals_fast(BL, FR);   % will pull from cached getters on client or you can pass BL/FR in
end

alpha = FRSTM_TH_PARAM.alpha;
M_oo = FRSTM_TH_PARAM.M_oo;
P_oo = FRSTM_TH_PARAM.P_oo;
T_oo = FRSTM_TH_PARAM.T_oo;
% M_th = 2.1;
m_dot = FRSTM_TH_PARAM.m_dot;
h_th = FRSTM_TH_PARAM.h_th;
T_th = FRSTM_TH_PARAM.T_th;
% PR_th = 120;

%% ---- compute -----

[obj_fn, ~, ~, ~] = FUN_objective_function(alpha, M_oo, P_oo, ...
    T_oo, m_dot, h_th, T_th, theta, 'n', 15);  % 15 dummy value to generate side fence

%% ---- Extract results ------
drag = obj_fn(1);
% l2d = obj_fn(2);
% pressureRatio = obj_fn(2);
pressureRecovery = obj_fn(2);
% intakeLength = obj_fn(5);
% intakeWidth = obj_fn(6);

% Weights
w_drag = 1;%0.01;
% w_l2d = 1;%0.01;
% w_PR  = 1;%0.5;
% w_len = 1;%0.01;
% w_wdt = 1;%0.01;
w_PRr = 1.0;

obj = [w_drag*drag, w_PRr*(-pressureRecovery)];
end
