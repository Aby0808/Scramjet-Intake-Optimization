function [ramp_coord1] = FUN_visc_blunt_correction(gamma, R, P_oo, T_oo,...
    ramp_param, ramp_coord, req_th_ht, mdot, flag, post)
% this function applies boundary layer and blunt edge correction to the
% inviscid ramp flow

mu_r = 1.79*10^-5;
ramp_param1 = ramp_param;

%% adjust ramp angles for boundary layer

for i=1:4
    
    if i==3
        l = sqrt((ramp_coord(4,1))^2 + (ramp_coord(4,2))^2);
    elseif i==4
        l = sqrt((ramp_coord(3,1)-ramp_coord(5,1))^2 + (ramp_coord(3,2)-ramp_coord(5,2))^2);
    else
        l = sqrt((ramp_coord(i,1)-ramp_coord(i+1,1))^2 + (ramp_coord(i,2)-ramp_coord(i+1,2))^2);
    end
    u = ramp_param(i,1)*sqrt(gamma*R*ramp_param(i,5));
    rho = ramp_param(i,4)/(R*ramp_param(i,5));
    mu = mu_r * ((ramp_param(i,5)/293.15)^1.5 * (403.15/(ramp_param(i,5)+110)));
    disp = ((rho*u/mu)^(-1/5) * (0.0156*l)^(4/5))/0.776;
    del = atand(disp/l);
    ramp_param1(i,2) = ramp_param1(i,2)-del;
end

ramp_param1(i+1,2) = ramp_param1(1,2)+ramp_param1(2,2)-ramp_param1(3,2)+ramp_param1(4,2);

rho_oo  = P_oo/(R*T_oo);
rho1 = ramp_param(1,4)/(R*ramp_param(1,5));
ssofd1 = 0.005*(rho_oo/rho1)/(1 + sqrt(rho_oo/rho1));

rho2  = ramp_param(2,4)/(R*ramp_param(2,5));
rho3 = ramp_param(3,4)/(R*ramp_param(3,5));
ssofd2 = 0.005*(rho2/rho3)/(1 + sqrt(rho2/rho3));

%% new ramp geoemtry

%% calculate cowl width

rho_cw = ramp_param(5,4)/(287*ramp_param(5,5));
V_cw = ramp_param(5,1)*sqrt(gamma*287*ramp_param(5,5));
wdt_cw = mdot/(rho_cw*V_cw*req_th_ht);

%% calculate ramp geometry

ramp_coord1 = zeros(5,2);
y_cap = 3;
th_diff = 10;

while th_diff>10^-5

    % first ramp
    m1 = tand(-ramp_param1(1,3));
    [x1,y1] = FUN_find_intersection(0,y_cap,0,-(ssofd1+ssofd2),-0.015,m1);
    ramp_coord1(1,:) = [x1,y1];

    % second ramp
    [x2,y2] = FUN_find_intersection(-(ssofd1+ssofd2),-0.015,tand(-(ramp_param1(1,2)+ramp_param1(2,3))), ...
        x1,y1,tand(-ramp_param1(1,2)));
    ramp_coord1(2,:) = [x2,y2];

    % third ramp
    [x3,y3] = FUN_find_intersection(-(ssofd1+ssofd2),ssofd2,tand(-(ramp_param1(1,2)+ramp_param1(2,2)-ramp_param1(3,3))),x2,y2, ...
        tand(-(ramp_param1(1,2)+ramp_param1(2,2))));
    ramp_coord1(3,:) = [x3,y3];

    % fourth ramp
    [x4,y4] = FUN_find_intersection(0,0.0,tand(ramp_param1(3,2)-ramp_param1(1,2)-ramp_param1(2,2)),x3,y3,...
        tand(ramp_param1(3,2)-ramp_param1(1,2)-ramp_param1(2,2)-ramp_param1(4,3)));
    ramp_coord1(4,:) = [x4,y4];

    % fifth ramp
    [x5,y5] = FUN_find_intersection(x3,y3,tand(ramp_param1(3,2)-ramp_param1(1,2)-ramp_param1(2,2)-ramp_param1(4,2)),x4,y4,...
        tand(ramp_param1(3,2)-ramp_param1(1,2)-ramp_param1(2,2)-ramp_param1(4,2)+ramp_param1(5,3)));
    ramp_coord1(5,:) = [x5,y5];

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
ln_in = ramp_coord1(5,1) - ramp_coord1(1,1);

%% plotting

if post == 'y'
    % figure
    plot([x1,x2,x3,x5],[y1,y2,y3,y5],Color='g')
    hold on
    plot([0,x4],[0.0,y4], Color='g')
    plot([x1,-(ssofd1+ssofd2),x2],[y1,-0.015,y2],Color='r')
    plot([-(ssofd1+ssofd2),x3,x4,x5],[0.005,y3,y4,y5],Color='r')
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