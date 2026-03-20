function [L,D,coord,int_wdt,ln_in] = FUN_compute_aero_forces_wr_sr(gamma, alpha,...
    M_oo, P_oo, ramp_coord, ramp_prop, cwl_wdt, post)

% this function computes and returns the drag and lift forces on the intake

%% intake width and length

% y1=ramp_coord(1,2)-tand(ramp_prop(1,2))*(ramp_coord(3,1)-ramp_coord(1,1));
% y2=ramp_coord(3,2);

sr_wdt=sqrt((ramp_coord(3,1)-ramp_coord(2,1))^2)*tand(ramp_prop(2,3));

% sr_wdt = sqrt((y2-y1)^2);

if (2*sr_wdt+cwl_wdt) > (1.5*cwl_wdt)
    int_wdt = 2*sr_wdt+cwl_wdt;
else
    int_wdt=1.5*cwl_wdt;
end
ln_in = ramp_coord(end-1,1) - ramp_coord(1,1);

%% first ramp using variable wedge waverider

[L1, D1,coord,ramp_ar1] = FUN_var_wdg_waverider_aero(gamma,alpha, P_oo, ramp_prop, ramp_coord, int_wdt, cwl_wdt);

%% second ramp

inc_l = ramp_prop(1,2)+ramp_prop(2,2);  % local inclination angle
ramp_ar2 = cwl_wdt*sqrt((ramp_coord(2,1)-ramp_coord(3,1))^2 ...
    + (ramp_coord(2,2)-ramp_coord(3,2))^2);
length2 = sqrt((ramp_coord(2,1)-ramp_coord(3,1))^2 + (ramp_coord(2,2)-ramp_coord(3,2))^2);
force = (ramp_prop(2,4) - P_oo)*ramp_ar2;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(2,1), ramp_prop(2,4), ...
    ramp_prop(2,5), 400, length2, ramp_ar2);
L2 = force*sind(90-inc_l) - force_v*sind(inc_l);
D2 = force*cosd(90-inc_l) + force_v*cosd(inc_l);

%% third ramp

% upper ramp
inc_l = ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,2);
length3 = sqrt((ramp_coord(3,1)-ramp_coord(4,1))^2 + (ramp_coord(3,2)-ramp_coord(4,2))^2);
ramp_ar3 = cwl_wdt * length3;
force = (ramp_prop(3,4) - P_oo)*ramp_ar3;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(3,1), ramp_prop(3,4), ...
    ramp_prop(3,5), 600, length3, ramp_ar3);
L3u = force*sind(90-inc_l) - force_v*sind(inc_l);
D3u = force*cosd(90-inc_l) + force_v*cosd(inc_l);

% lower ramp
length3 = sqrt((ramp_coord(5,1))^2 + (ramp_coord(5,2))^2);
ramp_ar3 = cwl_wdt * length3;
force = (ramp_prop(3,4) - ramp_prop(2,4)*FUN_oblique_shock_prop_ratio(gamma,ramp_prop(3,1)...
    ,ramp_prop(3,2),'PR'))*ramp_ar3;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(3,1), ramp_prop(3,4), ...
    ramp_prop(3,5), 600, length3, ramp_ar3);
L3l = -force*sind(90-inc_l) - force_v*sind(inc_l);
D3l = -force*cosd(90-inc_l) + force_v*cosd(inc_l);
L3 = L3u + L3l;
D3 = D3u + D3l;

%% fourth (cowl) ramp

% upper ramp
inc_l = ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,2)-ramp_prop(4,2);
length4 = sqrt((ramp_coord(4,1)-ramp_coord(6,1))^2 + (ramp_coord(4,2)-ramp_coord(6,2))^2);
ramp_ar4 = cwl_wdt * length4;
force = (ramp_prop(4,4) - P_oo)*ramp_ar4;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(4,1), ramp_prop(4,4), ...
    ramp_prop(4,5), 600, length4, ramp_ar4);
L4u = force*sind(90-inc_l) - force_v*sind(inc_l);
D4u = force*cosd(90-inc_l) + force_v*cosd(inc_l);

% lower ramp
length4 = sqrt((ramp_coord(5,1)-ramp_coord(7,1))^2 + (ramp_coord(5,2)-ramp_coord(7,2))^2);
ramp_ar4 = cwl_wdt * length4;
force = (ramp_prop(4,4) - ramp_prop(3,4)*FUN_oblique_shock_prop_ratio(gamma,ramp_prop(2,1)...
    ,ramp_prop(3,2),'PR'))*ramp_ar4;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(4,1), ramp_prop(4,4), ...
    ramp_prop(4,5), 600, length4, ramp_ar4);
L4l = -force*sind(90-inc_l) - force_v*sind(inc_l);
D4l = -force*cosd(90-inc_l) + force_v*cosd(inc_l);
L4 = L4u + L4l;
D4 = D4u + D4l;

%% total

L1=L1*((ramp_ar1-ramp_ar2)/ramp_ar1);
D1=D1*((ramp_ar1-ramp_ar2)/ramp_ar1);  % correcting for area of forebody compression surface covered by second ramp

L = L1+L2+L3+L4;
D = D1+D2+D3+D4;

%% drag from blunt edges and side ramps

M2 = sqrt((1 + (0.5*(gamma-1))*M_oo^2)/(gamma*M_oo^2 - 0.5*(gamma-1))); % Mach number behind the normal shock
P_stag = P_oo*(1 + (2*gamma*(M_oo^2 -1)/(gamma+1)))*(1 + 0.5*(gamma-1)*M2^2)^(gamma/(gamma-1)); % stag pressure behind the normal shock

sr_nor_force = (ramp_prop(2,4)-P_oo)*(0.5*sqrt((ramp_coord(2,1)-ramp_coord(3,1))^2)*(int_wdt-cwl_wdt)); %total normal force on side ramp
sr_ax_force = (ramp_prop(2,4)-P_oo)*(0.5*sqrt((ramp_coord(2,2)-ramp_coord(3,2))^2)*(int_wdt-cwl_wdt)); %total axial force on side ramp

rot_mat = [cosd(alpha) sind(alpha);-sind(alpha) cosd(alpha)]; % rotation matrix for wind to body axis
fwa = rot_mat\[sr_nor_force; sr_ax_force]; % force in wind axis

sd_rp_ar = length2*sqrt(2)*(int_wdt-cwl_wdt);  % side ramp total area
Dv_sr = FUN_compute_viscous_drag_force(gamma, ramp_prop(2,1), ramp_prop(2,4),...
    ramp_prop(2,5), 600, length2, sd_rp_ar); % viscous force on side ramp

L_sr = fwa(1);          %net lift on side ramp
D_sr = fwa(2) + Dv_sr;  %net drag on side ramp

D_blunt = (P_stag - P_oo) * (int_wdt + cwl_wdt)*0.01;  % for 0.01m dia blunt on forebody leading edge and cowl edge

D = D + D_blunt + D_sr;
L = L + L_sr;

%% plotting step for postprocessing
if post == 'y'
    figure
    hold on
    view(3)
    grid on
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title('3D Intake Geometry')

    % Plot already plotted ramp 1 (points 1 and 2 are part of it)
    scatter3(coord(:,1), coord(:,2), coord(:,3), 'filled');

    % Compute intake width (int_wdt is top surface width, assumed symmetric about z=0)
    half_wdt = cwl_wdt / 2;

    % Define helper to plot planar surface
    plotRamp = @(p1, p2) fill3([p1(1), p2(1), p2(1), p1(1)], ...
                               [p1(2), p2(2), p2(2), p1(2)], ...
                               [-half_wdt, -half_wdt, half_wdt, half_wdt], ...
                               [0.8 0.8 1], 'FaceAlpha', 0.7, 'EdgeColor', 'k');

    % Plot ramps: 2→3, 3→4, 4→6 (upper ramps)
    for i = 2:3
        plotRamp(ramp_coord(i,:), ramp_coord(i+1,:));
    end
    plotRamp(ramp_coord(4,:), ramp_coord(6,:));

    plotRamp(ramp_coord(5,:), [0,0]);
    % Plot cowl: 5→7 (lower ramp)
    plotRamp(ramp_coord(5,:), ramp_coord(7,:));

    % Optionally: mark key points
    % for i = 1:7
    %     text(ramp_coord(i,1), ramp_coord(i,2), 0, num2str(i), ...
    %         'FontSize', 10, 'Color', 'r', 'HorizontalAlignment', 'center');
    % end
    axis equal
    hold off

end

end