function [L,D,coord,int_wdt,ln_in] = FUN_compute_aero_forces_wr_sf(alpha,...
    M_oo, P_oo, T_oo, ramp_coord, ramp_prop, cwl_wdt, post)

% this function computes and returns the drag and lift forces on the intake

%% first ramp using variable wedge waverider

[L1, D1,coord,ramp_ar1,int_wdt] = FUN_var_wdg_waverider_aero(alpha, P_oo, ramp_prop, ramp_coord, cwl_wdt);

%% second ramp

inc_l = ramp_prop(1,2)+ramp_prop(2,2);  % local inclination angle
ramp_ar2 = cwl_wdt*sqrt((ramp_coord(2,1)-ramp_coord(3,1))^2 ...
    + (ramp_coord(2,2)-ramp_coord(3,2))^2);
length2 = sqrt((ramp_coord(2,1)-ramp_coord(3,1))^2 + (ramp_coord(2,2)-ramp_coord(3,2))^2);
force = (ramp_prop(2,4) - P_oo)*ramp_ar2;
force_v = FUN_compute_viscous_drag_force(ramp_prop(2,1), ramp_prop(2,4), ...
    ramp_prop(2,5), length2, ramp_ar2);
L2 = force*cosd(inc_l+alpha) - force_v*sind(inc_l+alpha);
D2 = force*sind(inc_l+alpha) + force_v*cosd(inc_l+alpha);

%% third ramp

% upper ramp
inc_l = ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,2);
length3 = sqrt((ramp_coord(3,1)-ramp_coord(4,1))^2 + (ramp_coord(3,2)-ramp_coord(4,2))^2);
ramp_ar3 = cwl_wdt * length3;
force = (ramp_prop(3,4) - P_oo)*ramp_ar3;
force_v = FUN_compute_viscous_drag_force(ramp_prop(3,1), ramp_prop(3,4), ...
    ramp_prop(3,5), length3, ramp_ar3);
L3u = force*cosd(inc_l+alpha) - force_v*sind(inc_l+alpha);
D3u = force*sind(inc_l+alpha) + force_v*cosd(inc_l+alpha);

% lower ramp
length3 = sqrt((ramp_coord(5,1))^2 + (ramp_coord(5,2))^2);
ramp_ar3 = cwl_wdt * length3;
% force = (ramp_prop(3,4) - ramp_prop(2,4)*FUN_oblique_shock_prop_ratio(gamma,ramp_prop(3,1)...
%     ,ramp_prop(3,2),'PR'))*ramp_ar3;
force_v = FUN_compute_viscous_drag_force(ramp_prop(3,1), ramp_prop(3,4), ...
    ramp_prop(3,5), length3, ramp_ar3);
L3l = -force*cosd(inc_l+alpha) - force_v*sind(inc_l+alpha);
D3l = -force*sind(inc_l+alpha) + force_v*cosd(inc_l);
L3 = L3u + L3l;
D3 = D3u + D3l;

%% fourth (cowl) ramp

% upper ramp
inc_l = ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,2)-ramp_prop(4,2);
length4 = sqrt((ramp_coord(4,1)-ramp_coord(6,1))^2 + (ramp_coord(4,2)-ramp_coord(6,2))^2);
ramp_ar4 = cwl_wdt * length4;
force = (ramp_prop(4,4) - P_oo)*ramp_ar4;
force_v = FUN_compute_viscous_drag_force(ramp_prop(4,1), ramp_prop(4,4), ...
    ramp_prop(4,5), length4, ramp_ar4);
L4u = force*cosd(inc_l+alpha) - force_v*sind(inc_l+alpha);
D4u = force*sind(inc_l+alpha) + force_v*cosd(inc_l+alpha);

% lower ramp
length4 = sqrt((ramp_coord(5,1)-ramp_coord(7,1))^2 + (ramp_coord(5,2)-ramp_coord(7,2))^2);
ramp_ar4 = cwl_wdt * length4;
% force = (ramp_prop(4,4) - ramp_prop(3,4)*FUN_oblique_shock_prop_ratio(gamma,ramp_prop(2,1)...
    % ,ramp_prop(3,2),'PR'))*ramp_ar4;
force_v = FUN_compute_viscous_drag_force(ramp_prop(4,1), ramp_prop(4,4), ...
    ramp_prop(4,5), length4, ramp_ar4);
L4l = -force*cosd(inc_l+alpha) - force_v*sind(inc_l+alpha);
D4l = -force*sind(inc_l+alpha) + force_v*cosd(inc_l+alpha);
L4 = L4u + L4l;
D4 = D4u + D4l;

%% total

% L1=L1*((ramp_ar1-ramp_ar2)/ramp_ar1);
% D1=D1*((ramp_ar1-ramp_ar2)/ramp_ar1);  % correcting for area of forebody compression surface covered by second ramp

L = L1+L2+L3+L4;
D = D1+D2+D3+D4;

%% drag from blunt edges and side fences

% M2 = sqrt((1 + (0.5*(gamma-1))*M_oo^2)/(gamma*M_oo^2 - 0.5*(gamma-1))); % Mach number behind the normal shock
[ratio1, M2] = FUN_normal_shock_calculator_NASA9(M_oo, P_oo, T_oo);
% P_stag = P_oo*(1 + (2*gamma*(M_oo^2 -1)/(gamma+1)))*(1 + 0.5*(gamma-1)*M2^2)^(gamma/(gamma-1)); % stag pressure behind the normal shock
u = M2*sqrt(FUN_get_prop_NASA9(T_oo*ratio1(2),'gamma')*287*T_oo*ratio1(2));
[P_stag, ~] = FUN_get_stagnation_properties(T_oo*ratio1(2),P_oo*ratio1(1),u);

sd_fn = sqrt(ramp_coord(2,1)^2 + ramp_coord(2,2)^2); % ramp 2 side fence length
sd_fn_ar = 0.5*sqrt((ramp_coord(2,1)*ramp_coord(3,2) + ramp_coord(3,1)*(-ramp_coord(2,2)))^2); % side fence area using Heron's formula

Dv_sf = FUN_compute_viscous_drag_force(ramp_prop(2,1), ramp_prop(2,4), ramp_prop(2,5), sd_fn, sd_fn_ar);

D_sf = 2*(P_stag - P_oo)*sd_fn*0.01*cosd(ramp_prop(1,2)+ramp_prop(2,3)); % pressure drag on side fence
D_sf = D_sf + Dv_sf;
L_sf = 2*(P_stag - P_oo)*sd_fn*0.01*sind(ramp_prop(1,2)+ramp_prop(2,3)); %lift from side fence

D_blunt = (P_stag - P_oo) * (int_wdt + cwl_wdt)*0.01;  % for 0.01m dia blunt on forebody leading edge and cowl edge

D = D + D_blunt + D_sf;
L = L + L_sf;

%% intake width and length

% int_wdt = 1.5*cwl_wdt;
ln_in = ramp_coord(7,1) - coord(1,1);

%% plotting step for postprocessing
if post == 'y'
    figure
    hold on
    view(3)
    grid on
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title('3D Intake Geometry')

    % Re-orient first ramp point cloud: swap y and z
    scatter3(coord(:,1), coord(:,3), coord(:,2), 6, 'b', 'filled');

    % Compute intake width (int_wdt is top surface width, assumed symmetric about z=0)
    half_wdt = cwl_wdt / 2;

    % Define helper to plot planar ramp with axis swap
    plotRamp = @(p1, p2) fill3([p1(1), p2(1), p2(1), p1(1)], ...
                               [-half_wdt, -half_wdt, half_wdt, half_wdt], ...
                               [p1(2), p2(2), p2(2), p1(2)], ...
                               [0.8 0.8 1], 'FaceAlpha', 0.7, 'EdgeColor', 'k');

    % ramp1 waverider downstream
    zero_ind = find(coord(:,1)==0);
    if isempty(zero_ind)
        zero_ind = size(coord,1);
    end
    fill3([coord(zero_ind(1)-1,1),ramp_coord(2,1),ramp_coord(2,1),coord(zero_ind(1)-1,1)], ...
          [int_wdt/2,half_wdt,-half_wdt,-int_wdt/2], ...
          [coord(zero_ind(1)-1,2),ramp_coord(2,2),ramp_coord(2,2),coord(zero_ind(1)-1,2)], ...
          [0.8 0.8 1], 'FaceAlpha', 0.7, 'EdgeColor', 'k')

    % Plot ramps: 2→3, 3→4, 4→6 (upper ramps)
    for i = 2:3
        plotRamp(ramp_coord(i,:), ramp_coord(i+1,:));
    end
    plotRamp(ramp_coord(4,:), ramp_coord(6,:));

    % Plot cowl ramps: 5→7 and 5→(0,0)
    plotRamp(ramp_coord(5,:), [0, 0]);
    plotRamp(ramp_coord(5,:), ramp_coord(7,:));

    % Optional point labels
    % for i = 1:7
    %     text(ramp_coord(i,1), ramp_coord(i,2), 0, num2str(i), ...
    %         'FontSize', 10, 'Color', 'r', 'HorizontalAlignment', 'center');
    % end

    axis equal
    hold off
end

end