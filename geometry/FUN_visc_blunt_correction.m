function [ramp_param1, ramp_coord1, Tmdot, Mmdot, wdt_cwl] = FUN_visc_blunt_correction(M_oo,...
    ramp_param, ramp_coord, req_th_ht, mdot, post)

global BL_SHAPE_PARAM FRSTM_TH_PARAM

if isempty(BL_SHAPE_PARAM) || isempty(FRSTM_TH_PARAM)
    paths = setup_project_paths();
    BL = table2array(readtable(fullfile(paths.dataDir, 'bl_shap_param.dat')));
    FR = FUN_set_freestream_throat_params();
    FUN_setup_globals_fast(BL, FR);
end

P_oo = FRSTM_TH_PARAM.P_oo;
T_oo = FRSTM_TH_PARAM.T_oo;

R = 0.005; %blunt radius

% this function applies boundary layer and blunt edge correction to the inviscid ramp flow

mu_r = 1.79*10^-5;  %reference viscosity
get_mu = @(T) mu_r * ((T/293.15)^1.5 * (403.15/(T+110))); % dynamics viscosity using Sutherland's

ddel1 = @(r,u,m,x,h) h*(0.0125 * (r*u/m)^(-1/4) * x *(5/4))^(4/5);  %    0.215 * (r*u/m)^(5/4);   % change in momentum thickness

ramp_param1 = ramp_param; % copy to plot both old and new ramps
del = zeros(1,4);

% mdot = mdot1;

disp3 = 0;
T3mid = ramp_param1(3,5);
T4mid = ramp_param1(4,5);
M3mid = ramp_param1(3,2);
M4mid = ramp_param1(4,2);
P3mid = ramp_param1(3,4);
P4mid = ramp_param1(4,4);

l2_norm = 100;
max_iter = 1000;
n = 1;

% old values at throat
Mold = ramp_param(5,1);
Pold = ramp_param(5,4);
Told = ramp_param(5,5);

while l2_norm > 10^-5  % outer loop until flow properties at throat do not change


    %% Calculate the deflection cause from boundary layer
    for i=1:4

        if i==3
            l = 0.15; %sqrt((ramp_coord(4,1))^2 + (ramp_coord(4,2))^2);
            gamma = FUN_get_prop_NASA9(T3mid,'gamma');
            u = M3mid*sqrt(gamma*287*T3mid);
            rho = P3mid/(287*T3mid);
            mu = get_mu(T3mid);
        elseif i==4
            % del3 = atand(disp/l);
            l = 0.15; %sqrt((ramp_coord(3,1)-ramp_coord(5,1))^2 + (ramp_coord(3,2)-ramp_coord(5,2))^2);
            gamma = FUN_get_prop_NASA9(T4mid,'gamma');
            u = M4mid*sqrt(gamma*287*T4mid);
            rho = P4mid/(287*T4mid);
            mu = get_mu(T4mid);
            disp3 = disp;
        else
            l = sqrt((ramp_coord(i,1)-ramp_coord(i+1,1))^2 + (ramp_coord(i,2)-ramp_coord(i+1,2))^2);
            gamma = FUN_get_prop_NASA9(ramp_param1(i,5),'gamma');
            u = ramp_param1(i,1)*sqrt(gamma*287*ramp_param1(i,5));
            rho = ramp_param1(i,4)/(287*ramp_param1(i,5));
            mu = get_mu(ramp_param1(i,5));
        end
        % gamma = FUN_get_prop_NASA9(ramp_param1(i,5),'gamma');
        % u = ramp_param1(i,1)*sqrt(gamma*287*ramp_param1(i,5));
        % rho = ramp_param1(i,4)/(287*ramp_param1(i,5));
        % mu = get_mu(ramp_param1(i,5));
        H = interp1(BL_SHAPE_PARAM(:,1),BL_SHAPE_PARAM(:,i+1),ramp_param1(i,1));
        disp = ddel1(rho,u,mu,l,H);
        del(i) = atand(disp/l);
        % ramp_param1(i,2) = ramp_param1(i,2)-del;
    end
    disp4 = disp;

    %% Calculate shock stand-off distances for le and cowl
    %  Billig, F. S., “Shock-Wave Shapes Around Spherical- and Cylindrical-Nosed Bodies,”
    
    del1 = R*0.386*exp(4.67/M_oo^2);
    Rc1 = R*1.386*exp(1.8/(M_oo-1)^0.75);
    ssofd1 = R+del1;
    beta1 = ramp_param(1,3);
    y1ht = sqrt((((R+del1)/(Rc1*cotd(beta1)))+1)^2 -1)*(Rc1/tand(beta1));

    del2 = R*0.386*exp(4.67/ramp_param(2,1)^2);
    Rc2 = R*1.386*exp(1.8/(ramp_param(2,1)-1)^0.75);
    ssofd2 = R+del2;
    beta2 = ramp_param(3,3);
    y2ht = sqrt((((R+del2)/(Rc2*cotd(beta2)))+1)^2 -1)*(Rc2/tand(beta2));

    %% calculate ramp geometry

    ramp_coord1 = zeros(7,2);
    y_cap = 3;
    th_diff = 10;

    while th_diff>10^-5

        % first ramp
        m1 = tand(-ramp_param1(1,3));                            %slope of first shock
        %find intercept where first shock intersects a horizontal line kept at
        %a capture height(first shock has to impinge on cowl(0,0), this is for inviscid flow without blunt)
        %the shock is assumed to have a hyperbolic shape at the blunt
        [x1,y1] = FUN_find_intersection(0,y_cap-y1ht,0,-ssofd2,-(0.005+y2ht),m1);
        [~, ~, Tbluntle, Mbluntle, ~] = FUN_blunt_body_shock(0.005, ramp_param1(1,3),...
            M_oo, P_oo, T_oo, 0, -0.005, sum(ramp_param1(1:2,1)), x1, y1, 'le');  % get the 
        % change in temperature and mach number
        ramp_coord1(1,:) = [x1,y_cap]; %[x1,y1+(2*(ssofd2+0.005)^2 - 0.005)];

        % second ramp
        % find intercept where second shock intersects the first ramp
        % [x2,y2] = FUN_find_intersection(-(ssofd2+0.005),-0.01,tand(-(ramp_param1(1,2)+ramp_param1(2,3))) ...
        %     ,x1,y1,tand(-ramp_param1(1,2)));
        [x2,y2] = FUN_find_intersection(-ssofd2,-(y2ht+0.005),tand(-(ramp_param1(1,2)+ramp_param1(2,3))) ...
            ,x1,y1,tand(-ramp_param1(1,2)+del(1)));
        ramp_coord1(2,:) = [x2,y2];

        % third ramp
        % top part of third ramp
        % find inetercept where second ramp intersects cowl shock
        % [x3,y3] = FUN_find_intersection(0,(2*(ssofd2+0.005) - 0.005),tand(-(ramp_param1(1,2)+ramp_param1(2,2) ...
        %     -ramp_param1(3,3))),x2,y2,tand(-(ramp_param1(1,2)+ramp_param1(2,2))));
        % [x3,y3] = FUN_find_intersection(-ssofd2,y2ht-0.005,tand(-(ramp_param1(1,2)+ramp_param1(2,2) ...
        %     -ramp_param1(3,3))),x2,y2,tand(-(ramp_param1(1,2)+ramp_param1(2,2)-del(2))));

        % since this point is close to the curved shock, using asymptotic
        % relations will not work
        [coord_inter, coord_shockcw, Tbluntcw, Mbluntcw, ~] = FUN_blunt_body_shock(0.005,...
            ramp_param1(3,3),ramp_param1(2,1),ramp_param1(2,4),ramp_param1(2,5),...
            0, -0.005, sum(ramp_param1(1:2,2)), x2, y2, 'cw');  % get the change in temperature and mach number
        x3=coord_inter(1); y3=coord_inter(2);
        ramp_coord1(3,:) = [x3,y3];

        % top part of third ramp
        % find point along the top part of third ramp 15cm away from (x3,y3)
        m = tand(ramp_param1(3,2)-ramp_param1(1,2)-ramp_param1(2,2)+del(3));
        x4 = x3 + 0.15*sqrt(1/(1+m^2));
        y4 = y3 + 0.15*m*sqrt(1/(1+m^2));
        ramp_coord1(4,:) = [x4,y4];

        % bottom part of third ramp
        % find intercept of the intersection of bottom part of third ramp and the shock impinging on (x4,y4)
        m1 = tand(ramp_param1(3,2)+ramp_param1(4,3)-ramp_param1(1,2)-ramp_param1(2,2));
        m2 = tand(ramp_param1(3,2)-ramp_param1(1,2)-ramp_param1(2,2)-del(3));
        [x5,y5] = FUN_find_intersection(x4,y4,m1,0,0,m2);
        ramp_coord1(5,:) = [x5,y5];

        % fourth ramp
        % top part of fourth ramp
        % find point along the top part of fourth ramp 15cm away from (x5,y5)
        m = tand(ramp_param1(3,2)+ramp_param1(4,2)-ramp_param1(1,2)-ramp_param1(2,2)+del(4));
        x6 = x4 + 0.15*sqrt(1/(1+m^2));
        y6 = y4 + 0.15*m*sqrt(1/(1+m^2));
        ramp_coord1(6,:) = [x6,y6];

        % bottom part of fourth ramp
        % find intercept of the intersection of bottom part of fourth ramp and the shock impinging on (x6,y6)
        m1 = tand(ramp_param1(3,2)+ramp_param1(4,2)+ramp_param1(5,3)-ramp_param1(1,2)-ramp_param1(2,2));
        m2 = tand(ramp_param1(3,2)+ramp_param1(4,2)-ramp_param1(1,2)-ramp_param1(2,2)-1*del(4));
        [x7,y7] = FUN_find_intersection(x6,y6,m1,x5,y5,m2);
        ramp_coord1(7,:) = [x7,y7];

        ht_th = y6-y7;
        th_diff = sqrt((ht_th-req_th_ht)^2);                % difference between calculated and desired cowl height
        y_cap = y_cap + 0.1*(req_th_ht-ht_th)/req_th_ht;    % change capture height based on the difference

    end

    %% Calculate change in flow properties due to skin friction and area change due to bl for ramps 3 and 4

    % Generalised one dimensional flow relations

    % calculate cowl width
    rho_cw = ramp_param1(5,4)/(287*ramp_param1(5,5));
    V_cw = ramp_param1(5,1)*sqrt(FUN_get_prop_NASA9(ramp_param1(5,5),'gamma')*287*ramp_param1(5,5));
    wdt_cwl = mdot/(rho_cw*V_cw*req_th_ht);

    %calculate properties before the fourth shock
    l3 = (0.15+sqrt(x5^2 + y5^2))/2;
    M31 = ramp_param1(3,1);
    P31 = ramp_param1(3,4);
    T31 = ramp_param1(3,5);
    rho3 = ramp_param1(3,4)/(287*ramp_param1(3,5));
    v3 = ramp_param1(3,1)*sqrt(FUN_get_prop_NASA9(ramp_param1(3,5),'gamma')*287*ramp_param1(3,5));
    A31 = mdot/(rho3*v3);                                           % Area behind the third shock
    A32 = A31*(wdt_cwl - 2*disp3)/wdt_cwl;                            % Area before the fourth shock
    D31 = 2*A31/(wdt_cwl + (A31/wdt_cwl));                            % hydraulic dia behind third shock
    D32 = 2*A32/((wdt_cwl - 2*disp3) + (A31/(wdt_cwl - 2*disp3)));    % hydraulic dia before fourth shock
    [M32,P32,T32,P032] = FUN_compute_generalized_1d_flow(M31,P31,T31,l3,D31,D32,A31,A32);
    M3mid = 0.5*(M31 + M32);
    P3mid = 0.5*(P31 + P32);
    T3mid = 0.5*(T31 + T32);

    % recalculate oblique shock solution for the new flow properties
    [result] = FUN_oblique_shock_calculator_NASA9(ramp_param1(4,2), M32, P32, T32, 'beta');
    beta4 = result(1);
    ratio1 = result(2:end-1);
    M4 = result(end);
    P4 = ratio1(1)*P32;
    T4 = ratio1(2)*T32;
    P04 = ratio1(3)*P032;
    ramp_param1(4,:) = [M4, ramp_param1(4,2), beta4, P4, T4, P04];


    %calculate properties before the fifth shock
    l4 = (0.15+sqrt((x7 -x5)^2 + (y7-y5)^2))/2;
    M41 = ramp_param1(4,1);
    P41 = ramp_param1(4,4);
    T41 = ramp_param1(4,5);
    rho4 = ramp_param1(4,4)/(287*ramp_param1(4,5));
    v4 = ramp_param1(4,1)*sqrt(FUN_get_prop_NASA9(ramp_param1(4,5),'gamma')*287*ramp_param1(4,5));
    A41 = mdot/(rho4*v4);                                           % Area behind the fourth shock
    A42 = A41*(wdt_cwl - 2*disp4)/wdt_cwl;                            % Area before the fifth shock
    D41 = 2*A41/(wdt_cwl + (A41/wdt_cwl));                            % hydraulic dia behind fourth shock
    D42 = 2*A42/((wdt_cwl - 2*disp4) + (A41/(wdt_cwl - 2*disp4)));    % hydraulic dia before fifth shock
    [M42,P42,T42,P042] = FUN_compute_generalized_1d_flow(M41,P41,T41,l4,D41,D42,A41,A42);
    M4mid = 0.5*(M41 + M42);
    P4mid = 0.5*(P41 + P42);
    T4mid = 0.5*(T41 + T42);

    % recalculate oblique shock solution for the new flow properties
    [result] = FUN_oblique_shock_calculator_NASA9(ramp_param1(5,2), M42, P42, T42, 'beta');
    beta5 = result(1);
    ratio1 = result(2:end-1);
    M5 = result(end);
    P5 = ratio1(1)*P42;
    T5 = ratio1(2)*T42;
    P05 = ratio1(3)*P042;
    ramp_param1(5,:) = [M5, ramp_param1(5,2), beta5, P5, T5, P05];

    %% calculate mass flow averaged throat properties
    gamma = FUN_get_prop_NASA9(T5,'gamma');
    Cp = FUN_get_prop_NASA9(T5,'cp');
    k = FUN_get_prop_NASA9(T5,'k');
    mu_r = 1.79*10^-5;

    Pr = mu_r*Cp/k;
    r = Pr^(1/3);  % recovery factor

    u = M5*sqrt(gamma*287*T5);
    hw = FUN_get_prop_NASA9(T5,'h') + r*0.5*u^2;
    Tw = FUN_get_T_from_h(hw);

    Tbl = @(Twl, Tool, y, del) (Tool - Twl)*(y/del)^(1/7) + Twl;
    ubl = @(uool, y, del) uool*(y/del)^(1/7);

    % Numerically integrate

    Tmdot = 0;
    umdot = 0;
    mdot1 = 0;
    blt = disp4/0.125;
    dy = blt/1000;
    for yl = 0:dy:blt
        Tmdot = Tmdot + P5*ubl(u,yl,blt)*dy*wdt_cwl/287;
        umdot = umdot + P5*ubl(u,yl,blt)^2 *dy*wdt_cwl/(287*Tbl(Tw,T5,yl,blt));
        mdot1 = mdot1 + P5*ubl(u,yl,blt)*dy*wdt_cwl/(287*Tbl(Tw,T5,yl,blt));
    end

    mdot1 = 2*(mdot1 + (P5*ubl(u,1,1)*(req_th_ht/2 - blt)*wdt_cwl)/(287*Tbl(Tw,T5,1,1)));
    Tmdot = 2*(Tmdot + (P5*ubl(u,1,1)*(req_th_ht/2 - blt)*wdt_cwl)/287)/mdot1;
    umdot = 2*(umdot + P5*ubl(u,1,1)^2 *(req_th_ht/2 - blt)*wdt_cwl/(287*Tbl(Tw,T5,1,1)))/mdot1;

    gmdot = FUN_get_prop_NASA9(Tmdot,'gamma');

    Mmdot = umdot/sqrt(gmdot*287*Tmdot);

    % add the effect of deteched shock on throat temperature and Mach number

    % Mmdot = Mmdot - (Mbluntcw*mdotbluntcw - Mbluntle*mdotbluntle)/(mdotbluntcw+mdotbluntle);
    % Tmdot = Tmdot + (Tbluntcw*mdotbluntcw + Tbluntle*mdotbluntle)/(mdotbluntcw+mdotbluntle);
    Mmdot = Mmdot - 0.02*Mbluntcw - 0.02*Mbluntle; % why these equations and what is the source of these?
    Tmdot = Tmdot + 0.02*Tbluntcw + 0.02*Tbluntle; % i don't know, they just work ¯\_(ツ)_/¯ (found emperically)


    %% calculate l2 norm and update the old flow property values

    l2_norm = sqrt((Mold-M5)^2 + (Pold-P5)^2 + (Told-T5)^2);
    Pold = P5;
    Told = T5;
    Mold = M5;

    if n>max_iter
        error('number of iterations exceeded max iteration limit of 1000') % don't run calculations till eternity
    end
    n=n+1;

end

% final cowl width
rho_cw = ramp_param1(5,4)/(287*Tmdot);
u_cw = Mmdot*sqrt(FUN_get_prop_NASA9(Tmdot,'gamma')*287*Tmdot);
wdt_cwl = mdot/(rho_cw*u_cw*0.054);


%% plotting

if post == 'y'
    % figure
    % plot([x1,x2,x3,x5],[y1,y2,y3,y5],Color='g')
    % hold on
    % plot([0,x4],[0.0,y4], Color='g')
    % plot([x1,-(ssofd1+ssofd2),x2],[y1,-0.015,y2],Color='r')
    % plot([-(ssofd1+ssofd2),x3,x4,x5],[0.005,y3,y4,y5],Color='r')
    % xlabel('x (m)')
    % ylabel('y (m)')
    % grid on
    % axis equal

    figure
    plot([x1,x2,x3,x4,x6],[y1,y2,y3,y4,y6],Color='g')
    hold on
    plot([0,x5,x7],[0,y5,y7], Color='g')
    plot([x1-ssofd1,-1*(ssofd2),x2],[y1-y1ht,-1*(y2ht+0.005),y2],Color='r')
    % plot([-ssofd2,x3],[1*(y2ht-0.005),y3],Color='r')
    % plot([x1,-1*(ssofd2),x2],[y1-(2*ssofd1),-2*ssofd2,y2],Color='r')
    % plot([-ssofd2,x3],[2*ssofd2,y3],Color='r')
    plot([x4,x5],[y4,y5],Color='r')
    plot([x6,x7],[y6,y7],Color='r')
    % plot(coord_shockle(1,:),coord_shockle(2,:),Color='r')
    plot(coord_shockcw(1,:),coord_shockcw(2,:),Color='r')
    xlabel('x (m)')
    ylabel('y (m)')
    grid on
    axis equal
end
end


