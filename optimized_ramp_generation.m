
% this script generates the ramps from the optimized ramp
projectPaths = setup_project_paths();

global FRSTM_TH_PARAM

if isempty(FRSTM_TH_PARAM)
    FRSTM_TH_PARAM = FUN_set_freestream_throat_params();
end

if ~exist('x', 'var')
    load(fullfile(projectPaths.resultsDir, 'optimization_results_w_int_constraints.mat'), 'x');
end

% Constants
alpha = FRSTM_TH_PARAM.alpha;  % angle of attack
M_oo = FRSTM_TH_PARAM.M_oo;    % freestream mach
P_oo = FRSTM_TH_PARAM.P_oo;    % freestream pressure
T_oo = FRSTM_TH_PARAM.T_oo;    % freestream temperature
% M_th = 2.1;
m_dot = FRSTM_TH_PARAM.m_dot;  % required mass flow rate
h_th = FRSTM_TH_PARAM.h_th;    % required throat height
T_th = FRSTM_TH_PARAM.T_th;    % approx desired pressure ratio

index = 79;  % replace the number with the index of the configuration needed
ramp_theta = x(index,:);


[obj_fn_val, ramp_coord, coord_ramp1,~] = FUN_objective_function(alpha, M_oo, P_oo, T_oo,...
    m_dot, h_th, T_th, ramp_theta, 'y',10);


figure
scatter3(coord_ramp1(:,1),coord_ramp1(:,2),coord_ramp1(:,3))
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
axis equal

