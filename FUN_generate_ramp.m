function [L, D, ramp_param, ramp_coord, coord_ramp1, Tth, Mth, wdt_in, wdt_cw, ln_in, SI] = ...
    FUN_generate_ramp(alpha, M_oo,P_oo, T_oo, ramp_param, req_th_ht, mdot, post, config_param)

% this function generates the intake geometry based on the ramp parameters
% it checks if the max throat height condition is satisfied before
% returning a complete intake

% cowl lip is origin
% ramp param(i,j) has j parameters(M, theta, beta, P, T, P0) for i ramps


%% calculate ramp geometry

ramp_coord = zeros(7,2);
y_cap = 3;      % some initial value
th_diff = 10;   % some initial value

while th_diff>10^-5  % iterate until difference between the throat height gets close tot he required one 

    % first ramp
    m1 = tand(-ramp_param(1,3));                            %slope of first shock
    %find intercept where first shock intersects a horizontal line kept at
    %a capture height(first shock has to impinge on cowl(0,0), this is for inviscid flow without blunt)
    [x1,y1] = FUN_find_intersection(0,y_cap,0,0,0,m1);
    ramp_coord(1,:) = [x1,y1];

    % second ramp
    % find intercept where second shock intersects the first ramp
    [x2,y2] = FUN_find_intersection(0,0,tand(-(ramp_param(1,2)+ramp_param(2,3))) ...
        ,x1,y1,tand(-ramp_param(1,2)));
    ramp_coord(2,:) = [x2,y2];

    % third ramp
    % top part of third ramp
    % find inetercept where second ramp intersects cowl shock
    [x3,y3] = FUN_find_intersection(0,0,tand(-(ramp_param(1,2)+ramp_param(2,2) ...
        -ramp_param(3,3))),x2,y2,tand(-(ramp_param(1,2)+ramp_param(2,2))));
    ramp_coord(3,:) = [x3,y3];

    % top part of third ramp
    % find point along the top part of third ramp 15cm away from (x3,y3)
    m = tand(ramp_param(3,2)-ramp_param(1,2)-ramp_param(2,2));
    x4 = x3 + 0.15*sqrt(1/(1+m^2));
    y4 = y3 + 0.15*m*sqrt(1/(1+m^2));
    ramp_coord(4,:) = [x4,y4];

    % bottom part of third ramp
    % find intercept of the intersection of bottom part of third ramp and the shock impinging on (x4,y4)
    m1 = tand(ramp_param(3,2)+ramp_param(4,3)-ramp_param(1,2)-ramp_param(2,2));
    m2 = tand(ramp_param(3,2)-ramp_param(1,2)-ramp_param(2,2));
    [x5,y5] = FUN_find_intersection(x4,y4,m1,0,0,m2);
    ramp_coord(5,:) = [x5,y5];

    % fourth ramp
    % top part of fourth ramp
    % find point along the top part of fourth ramp 15cm away from (x5,y5)
    m = tand(ramp_param(3,2)+ramp_param(4,2)-ramp_param(1,2)-ramp_param(2,2));
    x6 = x4 + 0.15*sqrt(1/(1+m^2));
    y6 = y4 + 0.15*m*sqrt(1/(1+m^2));
    ramp_coord(6,:) = [x6,y6];

    % bottom part of fourth ramp
    % find intercept of the intersection of bottom part of fourth ramp and the shock impinging on (x6,y6)
    m1 = tand(ramp_param(3,2)+ramp_param(4,2)+ramp_param(5,3)-ramp_param(1,2)-ramp_param(2,2));
    m2 = tand(ramp_param(3,2)+ramp_param(4,2)-ramp_param(1,2)-ramp_param(2,2));
    [x7,y7] = FUN_find_intersection(x6,y6,m1,x5,y5,m2);
    ramp_coord(7,:) = [x7,y7];

    ht_th = y6-y7;
    th_diff = sqrt((ht_th-req_th_ht)^2);                % difference between calculated and desired cowl height
    y_cap = y_cap + 0.1*(req_th_ht-ht_th)/req_th_ht;    % change capture height based on the difference

end


%% apply viscous corrections and blunting
% might have to implement this in some other part of the code
[ramp_param, ramp_coord1, Tth, Mth, wdt_cw] = FUN_visc_blunt_correction(M_oo, ramp_param, ...
    ramp_coord, req_th_ht, mdot, post);

%% startability index
M = ramp_param(1,2);
gamma = FUN_get_prop_NASA9(sum(ramp_param(2:5,5))/4,'gamma');
CRi = (1/M)*((2/(gamma+1))*(1 + 0.5*(gamma-1)*M^2))^(0.5*(gamma+1)/(gamma-1));
CRk = (((gamma+1)*M^2)/((gamma-1)*M^2 + 2))^0.5 * (((gamma+1)*M^2)/(2*gamma*M^2 ...
    - (gamma-1)))^(1/(gamma-1));
Acw = sqrt(x3^2 + y3^2);
Ath = sqrt((x7-x6)^2 + (y7-y6)^2);
CR = Acw/Ath;

SI = ((1/CR)-(1/CRi)/(1/CRk)-(1/CRi));

%% compute aerodynamic forces

if config_param>=0 && config_param<25
    [L,D, coord_ramp1, wdt_in, ln_in] = FUN_compute_aero_forces_wr_sf(alpha, M_oo, P_oo, T_oo, ramp_coord1,...
    ramp_param, wdt_cw, post);  % waverider+sidefence combination
elseif config_param>=25 && config_param<50
    [L,D, coord_ramp1, wdt_in, ln_in] = FUN_compute_aero_forces_wr_sr(gamma, alpha, M_oo, P_oo, ramp_coord,...
    ramp_param, wdt_cw, post);  % waverider+sideramp combination
elseif config_param>=50 && config_param<75
    delete(gcp('nocreate')); % shut down parellel pool
    fprintf('\nthis configuration is not yet implemented')
    error('\nconfiguration parameter %d out of range',config_param)
    % FUN_generate_ramp_sr_sf(gamma, alpha, M_oo, P_oo, ...
    % ramp_param, req_th_ht, mdot, post);  % sideramp+sidefence combination
elseif config_param>=75 && config_param<=100
    delete(gcp('nocreate')); % shut down parellel pool
    fprintf('\nthis configuration is not yet implemented')
    error('\nconfiguration parameter %d out of range',config_param)
    % FUN_generate_ramp_sr_sr(gamma, alpha, M_oo, P_oo, ...
    % ramp_param, req_th_ht, mdot, post);  % sideramp+sideramp combination
else
    delete(gcp('nocreate')); % shut down parellel pool
    fprintf('\nthis configuration is not yet implemented')
    error('\nconfiguration parameter %d out of range',config_param)
end

%% plotting step for postprocessing
if post == 'y'
    figure
    plot([x1,x2,x3,x4,x6],[y1,y2,y3,y4,y6],Color='k')
    hold on
    plot([0,x5,x7],[0,y5,y7], Color='k')
    plot([x1,0,x2],[y1,0,y2],Color='b')
    plot([0,x3],[0,y3],Color='b')
    plot([x4,x5],[y4,y5],Color='b')
    plot([x6,x7],[y6,y7],Color='b')
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
    axis equal
end

end