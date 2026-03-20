
% this script generates the ramps from the optimized ramp

global FRSTM_TH_PARAM

% Constants
alpha = FRSTM_TH_PARAM(1);  % angle of attack
M_oo = FRSTM_TH_PARAM(2);   % freestream mach
P_oo = FRSTM_TH_PARAM(3);   % freestream pressure
T_oo = FRSTM_TH_PARAM(4);   % freestream temperature
% M_th = 2.1;
m_dot = FRSTM_TH_PARAM(5);  % required mass flow rate
h_th = FRSTM_TH_PARAM(6);   % required throat height
T_th = FRSTM_TH_PARAM(7);   % approx desired pressure ratio

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