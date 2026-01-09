function [L,D,coord] = FUN_compute_aero_forces(gamma, M_oo, P_oo, T_oo, ramp_coord, ramp_prop, int_wdt, cwl_wdt, flag)

% this function computes and returns the drag and lift forces on the intake

%% first ramp using variable wedge waverider

[L1, D1,coord,ramp_ar1] = FUN_var_wdg_waverider_aero(gamma, M_oo, P_oo, ramp_prop, ramp_coord, int_wdt, cwl_wdt);
length1 = sqrt((ramp_coord(1,1)-ramp_coord(2,1))^2 + (ramp_coord(1,2)-ramp_coord(2,2))^2);
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(1,1), ramp_prop(1,4), ...
    ramp_prop(1,5), 300, length1, ramp_ar1);
L1 = L1 - force_v*sind(ramp_prop(1,2));
D1 = D1 + force_v*cosd(ramp_prop(1,2));

%% second ramp
% ramp geometry is created using machline cutting method

inc_l = ramp_prop(1,2)+ramp_prop(2,2);  % local inclination angle
ramp_ar2 = 0.5*(int_wdt + cwl_wdt)*sqrt((ramp_coord(2,1)-ramp_coord(3,1))^2 ...
    + (ramp_coord(2,2)-ramp_coord(3,2))^2);
length2 = sqrt((ramp_coord(2,1)-ramp_coord(3,1))^2 + (ramp_coord(2,2)-ramp_coord(3,2))^2);
force = (ramp_prop(2,4) - P_oo)*ramp_ar2;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(2,1), ramp_prop(2,4), ...
    ramp_prop(2,5), 400, length2, ramp_ar2);
L2 = force*sind(90-inc_l) - force_v*sind(inc_l);
D2 = force*cosd(90-inc_l) + force_v*cosd(inc_l);

%% third ramp

inc_l = ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,2)+ramp_prop(4,2);
ramp_ar3 = cwl_wdt * sqrt((ramp_coord(3,1)-ramp_coord(5,1))^2 + (ramp_coord(3,2)-ramp_coord(5,2))^2);
length3 = sqrt((ramp_coord(3,1)-ramp_coord(5,1))^2 + (ramp_coord(3,2)-ramp_coord(5,2))^2);
force = (ramp_prop(4,4) - P_oo)*ramp_ar3;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(4,1), ramp_prop(4,4), ...
    ramp_prop(4,5), 600, length3, ramp_ar3);
L3 = force*sind(90-inc_l) - force_v*sind(inc_l);
if inc_l>0
    D3 = force*cosd(90-inc_l) + force_v*cosd(inc_l);
else
    D3 = -force*cosd(90+inc_l) + force_v*cosd(inc_l);
end

%% fourth (cowl) ramp

inc_l = ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,2);
ramp_ar4 = cwl_wdt * sqrt(ramp_coord(4,1)^2 + ramp_coord(4,2)^2);
length4 = sqrt((ramp_coord(4,1))^2 + (ramp_coord(4,2))^2);
force = (ramp_prop(3,4) - P_oo)*ramp_ar4;
force_v = FUN_compute_viscous_drag_force(gamma, ramp_prop(3,1), ramp_prop(3,4), ...
    ramp_prop(3,5), 600, length4, ramp_ar4);
L4 = -force*sind(90-inc_l) - force_v*sind(inc_l);
if inc_l<0
    D4 = force*cosd(90+inc_l) + force_v*cosd(inc_l);
else
    D4 = -force*cosd(90-inc_l) + force_v*cosd(inc_l);
end

%% total

L = L1+L2+L3+L4;
D = D1+D2+D3+D4;

%% drag from blunt edges

M2 = sqrt((1 + (0.5*(gamma-1))*M_oo^2)/(gamma*M_oo^2 - 0.5*(gamma-1)));
P_stag = P_oo*(1 + (2*gamma*(M_oo^2 -1)/(gamma+1)))*(1 + 0.5*(gamma-1)*M2^2)^(gamma/(gamma-1));

if flag>=0 && flag <0.5
    cwl_fn = sqrt(ramp_coord(3,1)^2 + ramp_coord(3,2)^2); % cowl side fence
    D_sf = 2*(P_stag - P_oo)*cwl_fn*0.01*cosd(ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,3));
    L_sf = -2*(P_stag - P_oo)*cwl_fn*0.01*sind(ramp_prop(1,2)+ramp_prop(2,2)-ramp_prop(3,3));
elseif flag>=0.5 && flag <=1
    cwl_fn = sqrt(ramp_coord(2,1)^2 + ramp_coord(2,2)^2); % ramp 2 side fence
    D_sf = 2*(P_stag - P_oo)*cwl_fn*0.01*cosd(ramp_prop(1,2)+ramp_prop(2,3));
    L_sf = 2*(P_stag - P_oo)*cwl_fn*0.01*sind(ramp_prop(1,2)+ramp_prop(2,3));
end

D_blunt = (P_stag - P_oo) * (int_wdt + cwl_wdt)*0.01;  % for 0.01m blunt

D = D + D_blunt + D_sf;
L = L + L_sf;

end