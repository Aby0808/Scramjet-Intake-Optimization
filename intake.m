% MATLAB Script for a Five-Shock Scramjet Intake
clear

% Freestream conditions
M_oo = 6.5; % Fixed freestream Mach number
P_oo = 1071; % Fixed freestream pressure (Pa)
T_oo = 279; % Fixed freestream temperature (K)

% Throat conditions
M_th = 2; % Fixed throat Mach number
m_dot = 15; % Fixed mass flow rate (kg/s)
gamma = 1.4; % Ratio of specific heats
R = 287; % Specific gas constant (J/kg-K)

PR_th = 100;
TR_th = 1000/T_oo;
TPR = 0.4;

% Compute freestream properties
rho_inf = P_oo / (R * T_oo);
u_inf = M_oo * sqrt(gamma * R * T_oo);
P0_oo = P_oo * (1 + 0.5*(gamma-1)*M_oo^2)^(gamma/(gamma-1));

ramp_prop = zeros(5,6,100);
step_size = 1;

% Generate combinations of shock angles and filter invalid setups
c=1;
fprintf('\nS no.\tMach\tPR\tTR\tTPR');
% first ramp
beta1min = asind(1/M_oo); beta1max = beta1min + 20;
%for beta1 = beta1min+step_size:step_size:beta1max
beta1=10;
    M1 = FUN_oblique_shock(gamma, M_oo, beta1, 'beta2M2');
    theta1 = FUN_oblique_shock(gamma, M_oo, beta1, 'theta');
    P1 = P_oo * FUN_oblique_shock_prop_ratio(gamma, M_oo, theta1, 'PR');
    T1 = T_oo * FUN_oblique_shock_prop_ratio(gamma, M_oo, theta1, 'TR');
    P01 = P0_oo * FUN_oblique_shock_prop_ratio(gamma, M_oo, theta1,'TPR');
    ramp_prop(1,:,c) = [M1,theta1,beta1,P1,T1,P01];

    if P01/P0_oo > TPR %&& P1/P_oo < 2*PR_th
        % second ramp
        beta2min = asind(1/M1); beta2max = beta2min + 20;
        %for beta2 = beta2min+step_size:step_size:beta2max
        beta2 = 10;
            M2 = FUN_oblique_shock(gamma, M1, beta2, 'beta2M2');
            theta2 = FUN_oblique_shock(gamma, M1, beta2, 'theta');
            P2 = P1 * FUN_oblique_shock_prop_ratio(gamma, M1, theta2, 'PR');
            T2 = T1 * FUN_oblique_shock_prop_ratio(gamma, M1, theta2, 'TR');
            P02 = P01 * FUN_oblique_shock_prop_ratio(gamma, M1, theta2,'TPR');
            ramp_prop(2,:,c) = [M2,theta2,beta2,P2,T2,P02];

            if P02/P0_oo > TPR %&& P2/P_oo < 2*PR_th
                % third shock
                % turning back in because of straight cowl
                theta3min = (theta1+theta2)*0.5; theta3max = (theta1+theta2)*1.5;
                for theta3 = theta3min:step_size:theta3max
                    if theta3max < FUN_oblique_shock(gamma, M2, 0, 'theta_max')
                        M3 = FUN_oblique_shock(gamma, M2, theta3, 'theta2M2');
                        beta3 = FUN_oblique_shock(gamma, M2, theta3, 'beta');
                        P3 = P2 * FUN_oblique_shock_prop_ratio(gamma, M2, theta3, 'PR');
                        T3 = T2 * FUN_oblique_shock_prop_ratio(gamma, M2, theta3, 'TR');
                        P03 = P02 * FUN_oblique_shock_prop_ratio(gamma, M2, theta3,'TPR');
                        ramp_prop(3,:,c) = [M3,theta3,beta3,P3,T3,P03];

                        if P03/P0_oo > TPR %&& P3/P_oo < 2*PR_th
                            % fourth shock
                            beta4min = asind(1/M3); beta4max = beta4min + 30;
                            for beta4 = beta4min+step_size:step_size:beta4max
                                M4 = FUN_oblique_shock(gamma, M3, beta4, 'beta2M2');
                                theta4 = FUN_oblique_shock(gamma, M3, beta4, 'theta');
                                P4 = P3 * FUN_oblique_shock_prop_ratio(gamma, M3, theta4, 'PR');
                                T4 = T3 * FUN_oblique_shock_prop_ratio(gamma, M3, theta4, 'TR');
                                P04 = P03 * FUN_oblique_shock_prop_ratio(gamma, M3, theta4,'TPR');
                                ramp_prop(4,:,c) = [M4,theta4,beta4,P4,T4,P04];

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
                                        ramp_prop(5,:,c) = [M5,theta5,beta5,P5,T5,P05];

                                        % store and print filtered combination
                                        %if M5 > 1.5 && M5<2 && P05/P0_oo > TPR && P5/P_oo > 0.7*PR_th && P5/P_oo < 2*PR_th ...
                                                %&& T5/T_oo > 0.7*TR_th
                                            [L,D,ramp_coord,coord_ramp1,wdt_in,wdt_cwl,lng_in] = FUN_generate_ramp(gamma, M_oo,...
                                                P_oo, T_oo, ramp_prop, 0.055, 15);
                                            ob_fn = ((L/D)/D)*(P05/P0_oo)/(wdt_in*lng_in);
                                            fprintf('\n%i\t%f\t%f\t%f\t%f',c,M5,P5/P_oo,T5/T_oo,P05/P0_oo);
                                            c=c+1;

                                       % end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        %end
    end
%end


