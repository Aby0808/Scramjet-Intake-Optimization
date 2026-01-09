
% this script generates the ramps from the optimized ramp

% Constants
gamma = 1.4;
R = 287;
M_oo = 6.5;
P_oo = 1171;
T_oo = 279;
M_th = 2.5;
m_dot = 18.7;
PR_th = 100;

index = 74;
ramp_theta = x(index,:);


[obj_fn_val, ramp_coord, coord_ramp1] = FUN_objective_function(gamma, R, M_oo, P_oo, T_oo, M_th,...
    m_dot, PR_th, ramp_theta, x(index,end), 'y');


figure
scatter3(coord_ramp1(:,1),coord_ramp1(:,2),coord_ramp1(:,3))
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
axis equal