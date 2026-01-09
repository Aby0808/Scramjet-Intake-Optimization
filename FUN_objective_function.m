function [obj_fn, ramp_coord, coord_ramp1] = FUN_objective_function(gamma, R, M_oo, P_oo, T_oo, M_th, m_dot, ...
    PR_th, theta, flag, post)

% this function generates the objective function for teh optimization program

TR_th = 1030/T_oo;
TPR = 0.2;

penalty = [10^6, 10^6, 10^6, 10^2, 10^2, 10^2];

% Compute freestream properties
rho_inf = P_oo / (R * T_oo);
u_inf = M_oo * sqrt(gamma * R * T_oo);
P0_oo = P_oo * (1 + 0.5*(gamma-1)*M_oo^2)^(gamma/(gamma-1));

ramp_prop = zeros(5,6);

% Generate combinations of shock angles and filter invalid setups

% first ramp
theta1 = theta(1);
beta1=FUN_oblique_shock(gamma, M_oo, theta(1),'beta');
M1 = FUN_oblique_shock(gamma, M_oo, beta1, 'beta2M2');
P1 = P_oo * FUN_oblique_shock_prop_ratio(gamma, M_oo, theta1, 'PR');
T1 = T_oo * FUN_oblique_shock_prop_ratio(gamma, M_oo, theta1, 'TR');
P01 = P0_oo * FUN_oblique_shock_prop_ratio(gamma, M_oo, theta1,'TPR');
ramp_prop(1,:) = [M1,theta1,beta1,P1,T1,P01];

if P01/P0_oo > TPR && theta(2) < FUN_oblique_shock(gamma, M1, 0, 'theta_max')
    % second ramp
    theta2 = theta(2);
    beta2 = FUN_oblique_shock(gamma, M1, theta2, 'beta');
    M2 = FUN_oblique_shock(gamma, M1, beta2, 'beta2M2');
    P2 = P1 * FUN_oblique_shock_prop_ratio(gamma, M1, theta2, 'PR');
    T2 = T1 * FUN_oblique_shock_prop_ratio(gamma, M1, theta2, 'TR');
    P02 = P01 * FUN_oblique_shock_prop_ratio(gamma, M1, theta2,'TPR');
    ramp_prop(2,:) = [M2,theta2,beta2,P2,T2,P02];

    if P02/P0_oo > TPR && theta(3) < FUN_oblique_shock(gamma, M2, 0, 'theta_max')
        % third shock
        % turning back in because of straight cowl
        theta3 = theta(3);
        M3 = FUN_oblique_shock(gamma, M2, theta3, 'theta2M2');
        beta3 = FUN_oblique_shock(gamma, M2, theta3, 'beta');
        P3 = P2 * FUN_oblique_shock_prop_ratio(gamma, M2, theta3, 'PR');
        T3 = T2 * FUN_oblique_shock_prop_ratio(gamma, M2, theta3, 'TR');
        P03 = P02 * FUN_oblique_shock_prop_ratio(gamma, M2, theta3,'TPR');
        ramp_prop(3,:) = [M3,theta3,beta3,P3,T3,P03];

        if P03/P0_oo > TPR && theta(4) < FUN_oblique_shock(gamma, M3, 0, 'theta_max')
            % fourth shock
            theta4 = theta(4);
            beta4 = FUN_oblique_shock(gamma, M3, theta4, 'beta');
            M4 = FUN_oblique_shock(gamma, M3, beta4, 'beta2M2');
            theta4 = FUN_oblique_shock(gamma, M3, beta4, 'theta');
            P4 = P3 * FUN_oblique_shock_prop_ratio(gamma, M3, theta4, 'PR');
            T4 = T3 * FUN_oblique_shock_prop_ratio(gamma, M3, theta4, 'TR');
            P04 = P03 * FUN_oblique_shock_prop_ratio(gamma, M3, theta4,'TPR');
            ramp_prop(4,:) = [M4,theta4,beta4,P4,T4,P04];

            if P04/P0_oo > TPR %&& P4/P_oo < 2*PR_th
                % fifth shock
                % have to straighten the flow
                theta5 = theta1 + theta2 - theta3 + theta4;
                if theta5>0 && theta5<FUN_oblique_shock(gamma, M4, 0, 'theta_max')
                    beta5 = FUN_oblique_shock(gamma, M4, theta5, 'beta');
                    M5 = FUN_oblique_shock(gamma, M4, theta5, 'theta2M2');
                    P5 = P4 * FUN_oblique_shock_prop_ratio(gamma, M4, theta5, 'PR');
                    T5 = T4 * FUN_oblique_shock_prop_ratio(gamma, M4, theta5, 'TR');
                    P05 = P04 * FUN_oblique_shock_prop_ratio(gamma, M4, theta5,'TPR');
                    ramp_prop(5,:) = [M5,theta5,beta5,P5,T5,P05];

                    % store and print filtered combination
                    if P5/P_oo > 0.8*PR_th && P5/P_oo < 2*PR_th && P3/P2<3 ...
                            && T5/T_oo > 0.9*TR_th && T5/T_oo < 1.2*TR_th && P4/P3<1.7 && P5/P4<1.7
                        [L,D,ramp_coord,coord_ramp1,wdt_in,wdt_cwl,lng_in] = FUN_generate_ramp(gamma, M_oo,...
                            P_oo, T_oo, ramp_prop, 0.054, m_dot, flag, post);

                        if post=='y'
                            [~] = FUN_visc_blunt_correction(gamma, R, P_oo, T_oo, ramp_prop, ramp_coord, 0.064,m_dot,flag,post);
                        end

                        obj_fn = [D, L/D, P5/P_oo, P05/P0_oo, lng_in, wdt_in];
                        return
                    else
                        obj_fn = penalty;
                        ramp_coord = 0;
                        coord_ramp1 = 0;
                        return;
                    end
                else
                    obj_fn = penalty;
                    ramp_coord = 0;
                    coord_ramp1 = 0;
                    return;
                end
            else
                obj_fn = penalty;
                ramp_coord = 0;
                coord_ramp1 = 0;
                return;
            end
        else
            obj_fn = penalty;
            ramp_coord = 0;
            coord_ramp1 = 0;
            return;
        end
    else
        obj_fn = penalty;
        ramp_coord = 0;
        coord_ramp1 = 0;
        return;
    end
else
    obj_fn = penalty;
    ramp_coord = 0;
    coord_ramp1 = 0;
    return;
end

end