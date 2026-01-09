function [L, D, ramp_coord, coord_ramp1, wdt_in, wdt_cw, ln_in] = FUN_generate_ramp(gamma, M_oo, P_oo, T_oo,...
    ramp_param, req_th_ht, mdot, flag, post)

% this function generates the intake geometry based on the ramp parameters
% it checks if the max throat height condition is satisfied before
% returning a complete intake

% cowl lip is origin
% ramp param(i,j) has j parameters(M, theta, beta, P, T) for i ramps


%% calculate cowl width

rho_cw = ramp_param(5,4)/(287*ramp_param(5,5));
V_cw = ramp_param(5,1)*sqrt(gamma*287*ramp_param(5,5));
wdt_cw = mdot/(rho_cw*V_cw*req_th_ht);

%% calculate ramp geometry

ramp_coord = zeros(5,2);
y_cap = 3;
th_diff = 10;

while th_diff>10^-5

    % first ramp
    m1 = tand(-ramp_param(1,3));
    [x1,y1] = FUN_find_intersection(0,y_cap,0,0,0,m1);
    ramp_coord(1,:) = [x1,y1];

    % second ramp
    [x2,y2] = FUN_find_intersection(0,0,tand(-(ramp_param(1,2)+ramp_param(2,3))),x1,y1,tand(-ramp_param(1,2)));
    ramp_coord(2,:) = [x2,y2];

    % third ramp
    [x3,y3] = FUN_find_intersection(0,0,tand(-(ramp_param(1,2)+ramp_param(2,2)-ramp_param(3,3))),x2,y2, ...
        tand(-(ramp_param(1,2)+ramp_param(2,2))));
    ramp_coord(3,:) = [x3,y3];

    % fourth ramp
    [x4,y4] = FUN_find_intersection(0,0,tand(ramp_param(3,2)-ramp_param(1,2)-ramp_param(2,2)),x3,y3,...
        tand(ramp_param(3,2)-ramp_param(1,2)-ramp_param(2,2)-ramp_param(4,3)));
    ramp_coord(4,:) = [x4,y4];

    % fifth ramp
    [x5,y5] = FUN_find_intersection(x3,y3,tand(ramp_param(3,2)-ramp_param(1,2)-ramp_param(2,2)-ramp_param(4,2)),x4,y4,...
        tand(ramp_param(3,2)-ramp_param(1,2)-ramp_param(2,2)-ramp_param(4,2)+ramp_param(5,3)));
    ramp_coord(5,:) = [x5,y5];

    ht_th = y5-y4;
    th_diff = sqrt((ht_th-req_th_ht)^2);
    y_cap = y_cap + 0.1*(req_th_ht-ht_th)/req_th_ht;

end

%% intake width and length
mu2 = asind(1/ramp_param(2,1));
if flag>=0 && flag <0.5
    wdt_in = wdt_cw + 2*sqrt((x3-x2)^2 + (y3-y2)^2)*tand(mu2);
elseif flag <=1 && flag>=0.5
    wdt_in = 1.5*wdt_cw;
end
ln_in = ramp_coord(5,1) - ramp_coord(1,1);


%% apply viscous corrections and blunting


%% compute aerodynamic forces

% [L,D,coord_ramp1] = FUN_var_wdg_waverider_aero(gamma, M_oo, ramp_param, ramp_coord, wdt_in, wdt_cw);
[L,D, coord_ramp1] = FUN_compute_aero_forces(gamma, M_oo, P_oo, T_oo, ramp_coord, ramp_param,wdt_in, wdt_cw, flag);

% 

if post == 'y'
    figure
    plot([x1,x2,x3,x5],[y1,y2,y3,y5],Color='k')
    hold on
    plot([0,x4],[0,y4], Color='k')
    plot([x1,0,x2],[y1,0,y2],Color='b')
    plot([0,x3,x4,x5],[0,y3,y4,y5],Color='b')
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
    axis equal

    f_tip = wdt_in*(1 + 2*tand(asind(1/ramp_param(1,1))));

    fprintf('\nlength %d',ln_in)
    fprintf('\nwidth %d',wdt_in)
    fprintf('\ncowl %d',wdt_cw)
    fprintf('\nreuired flat portion length of tip %d\n',f_tip)
end

end