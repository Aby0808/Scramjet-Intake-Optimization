function [result] = FUN_conical_shock(gamma, M, angle, parameter)
% This function solves the Taylor-Maccoll equations to solve flow over a cone.
% Depending on 'parameter', it returns either the shock angle ('beta'),
% deflection angle ('theta'), or post-cone Mach number ('theta2M2').

switch parameter

    case 'beta'
        % Case to find shock angle beta
        beta = FUN_oblique_shock(gamma, M, angle, 'beta');

        beta_cone = beta;
        thetaf = 100;
        while sqrt((thetaf - angle)^2) >= 0.01
            M2 = FUN_oblique_shock(gamma, M, beta_cone, 'beta2M2');
            prop = zeros(100, 4);
            theta_cone = FUN_oblique_shock(gamma, M, beta_cone, 'theta');
            prop(1,1) = cosd(beta_cone - theta_cone) * M2;
            prop(1,2) = -sind(beta_cone - theta_cone) * M2;
            prop(1,3) = beta_cone;
            r = zeros(100,1);
            r(1) = 1;
            prop(1,4) = r(1);
            dtheta = -0.01;
            i = 1;
            while prop(i,2) <= 0
                temp = (prop(i,1) + prop(i,2)*cotd(prop(i,3))) / (prop(i,2)^2 - 1);
                dudtheta = prop(i,2) + 0.5*(gamma-1)*prop(i,1)*prop(i,2)*temp;
                dvdtheta = -prop(i,1) + (1 + 0.5*(gamma-1)*prop(i,2)^2)*temp;
                prop(i+1,1) = prop(i,1) + dudtheta*(dtheta*pi/180);
                prop(i+1,2) = prop(i,2) + dvdtheta*(dtheta*pi/180);
                prop(i+1,3) = prop(i,3) + dtheta;
                r(i+1) = r(i) + r(i)*prop(i,1)*(dtheta*pi/180)/prop(i,2);
                prop(i,4) = r(i+1);
                i = i + 1;
            end
            thetaf = prop(i,3);
            beta_cone = beta_cone + (angle - thetaf) / angle;
        end
        result = beta_cone;

    case 'theta'
        % Case to find final deflection angle theta
        beta_cone = angle;
        M2 = FUN_oblique_shock(gamma, M, beta_cone, 'beta2M2');
        prop = zeros(100, 4);
        theta_cone = FUN_oblique_shock(gamma, M, beta_cone, 'theta');
        prop(1,1) = cosd(beta_cone - theta_cone) * M2;
        prop(1,2) = -sind(beta_cone - theta_cone) * M2;
        prop(1,3) = beta_cone;
        r = zeros(100,1);
        r(1) = 1;
        prop(1,4) = r(1);
        dtheta = -0.001;
        i = 1;
        while prop(i,2) <= 0
            temp = (prop(i,1) + prop(i,2)*cotd(prop(i,3))) / (prop(i,2)^2 - 1);
            dudtheta = prop(i,2) + 0.5*(gamma-1)*prop(i,1)*prop(i,2)*temp;
            dvdtheta = -prop(i,1) + (1 + 0.5*(gamma-1)*prop(i,2)^2)*temp;
            prop(i+1,1) = prop(i,1) + dudtheta*(dtheta*pi/180);
            prop(i+1,2) = prop(i,2) + dvdtheta*(dtheta*pi/180);
            prop(i+1,3) = prop(i,3) + dtheta;
            r(i+1) = r(i) + r(i)*prop(i,1)*(dtheta*pi/180)/prop(i,2);
            prop(i,4) = r(i+1);
            i = i + 1;
        end
        thetaf = prop(i,3);
        result = thetaf;

    case 'theta2M2'
        % case returns final Mach number based on integration
        beta = FUN_oblique_shock(gamma, M, angle, 'beta');

        beta_cone = beta;
        thetaf = 100;
        while sqrt((thetaf - angle)^2) >= 0.01
            M2 = FUN_oblique_shock(gamma, M, beta_cone, 'beta2M2');
            prop = zeros(100, 4);
            theta_cone = FUN_oblique_shock(gamma, M, beta_cone, 'theta');
            prop(1,1) = cosd(beta_cone - theta_cone) * M2;
            prop(1,2) = -sind(beta_cone - theta_cone) * M2;
            prop(1,3) = beta_cone;
            r = zeros(100,1);
            r(1) = 1;
            prop(1,4) = r(1);
            dtheta = -0.01;
            i = 1;
            while prop(i,2) <= 0
                temp = (prop(i,1) + prop(i,2)*cotd(prop(i,3))) / (prop(i,2)^2 - 1);
                dudtheta = prop(i,2) + 0.5*(gamma-1)*prop(i,1)*prop(i,2)*temp;
                dvdtheta = -prop(i,1) + (1 + 0.5*(gamma-1)*prop(i,2)^2)*temp;
                prop(i+1,1) = prop(i,1) + dudtheta*(dtheta*pi/180);
                prop(i+1,2) = prop(i,2) + dvdtheta*(dtheta*pi/180);
                prop(i+1,3) = prop(i,3) + dtheta;
                r(i+1) = r(i) + r(i)*prop(i,1)*(dtheta*pi/180)/prop(i,2);
                prop(i,4) = r(i+1);
                i = i + 1;
            end
            thetaf = prop(i,3);
            beta_cone = beta_cone + (angle - thetaf) / angle;
        end
        % Final velocity components
        u_final = prop(i,1);
        v_final = prop(i,2);

        % Final Mach number from u and v
        M_final = sqrt(u_final^2 + v_final^2);

        result = M_final;

    otherwise
        error('Invalid parameter requested: %s', parameter);

end

end
