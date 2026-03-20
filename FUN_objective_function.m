function [obj_fn, ramp_coord, coord_ramp1, postprocess] = FUN_objective_function(alpha, M_oo, P_oo, ...
    T_oo, m_dot, h_th, T_t, theta, post, config_param)

% this function generates the objective function for the optimization program
% currently it generates two external and two internal ramps (two external and three internal shocks)
% first the ramps are generated for 2d flowfield
% the ramp combinations are filtered out and sent for viscous and blunt corrections, 3d ramp genration 
% and computing of aero forces
% the cart. coord. system is alinged with the body reference frame and not the wind refernce frame
% incoming flow is asumed to be inclined at angle alpha to the x axis of our coord. system
% first ramp is assumed to be alinged with the x axis and the subsequent ramps are alinged with respect
% to the previous ramp

%% initializing required parameters
% TR_th = T_t/T_oo;
gamma = FUN_get_prop_NASA9(T_oo,'gamma');

% Compute freestream properties
u_oo = M_oo*sqrt(FUN_get_prop_NASA9(T_oo,'gamma')*287*T_oo);
[P0_oo, ~] = FUN_get_stagnation_properties(T_oo,P_oo,u_oo);
% P0_oo = P_oo * (1 + 0.5*(gamma-1)*M_oo^2)^(gamma/(gamma-1));  % for some reason the 
% stagnation properties are not coming properly

penalty = [10^5, -10^5];  % penalty to be returned for undesired flowfield
postprocess = [0, 0, 0, 0, 0, 0];

ramp_prop = zeros(5,6);

try
    %% ramp generation and filtering

    % Iterate on the received ramp angles
    % Generate combinations of shock angles and filter invalid setups

    % net deflection from first two ramp should be greater than defplection of third ramp
    if theta(1)+theta(2)-theta(3)>0 && theta(1)+theta(2)-theta(3)-theta(4)>0 ...
            && theta(1) < FUN_oblique_shock(gamma, M_oo, 0, 'theta_max')

        % first ramp
        theta1 = theta(1);   % preserve shock strength and align beta back to coord system
        [result] = FUN_oblique_shock_calculator_NASA9(theta1+alpha, M_oo, P_oo, T_oo, 'beta');
        beta1 = result(1) - alpha;
        ratio1 = result(2:end-1);
        M1 = result(end);
        P1 = ratio1(1)*P_oo;
        T1 = ratio1(2)*T_oo;
        P01 = ratio1(3)*P0_oo;
        % P01 = P1 * (1 + 0.5*(FUN_get_prop_NASA9(T1,'gamma')-1)*M1^2)^...
        %    (FUN_get_prop_NASA9(T1,'gamma')/(FUN_get_prop_NASA9(T1,'gamma')-1));
        ramp_prop(1,:) = [M1,theta1,beta1,P1,T1,P01];

        % filtering based on max rmap angle
        gamma = FUN_get_prop_NASA9(T1,'gamma');
        if theta(2) < FUN_oblique_shock(gamma, M1, 0, 'theta_max')
            % second ramp
            theta2 = theta(2);
            [result] = FUN_oblique_shock_calculator_NASA9(theta2, M1, P1, T1, 'beta');
            beta2 = result(1);
            ratio1 = result(2:end-1);
            M2 = result(end);
            P2 = ratio1(1)*P1;
            T2 = ratio1(2)*T1;
            P02 = ratio1(3)*P01;
            % P02 = P2 * (1 + 0.5*(FUN_get_prop_NASA9(T2,'gamma')-1)*M2^2)^...
            %     (FUN_get_prop_NASA9(T2,'gamma')/(FUN_get_prop_NASA9(T2,'gamma')-1));
            ramp_prop(2,:) = [M2,theta2,beta2,P2,T2,P02];

            gamma = FUN_get_prop_NASA9(T2,'gamma');
            if theta(3) < FUN_oblique_shock(gamma, M2, 0, 'theta_max')
                % third shock
                theta3 = theta(3);
                [result] = FUN_oblique_shock_calculator_NASA9(theta3, M2, P2, T2, 'beta');
                beta3 = result(1);
                ratio1 = result(2:end-1);
                M3 = result(end);
                P3 = ratio1(1)*P2;
                T3 = ratio1(2)*T2;
                P03 = ratio1(3)*P02;
                % P03 = P3 * (1 + 0.5*(FUN_get_prop_NASA9(T3,'gamma')-1)*M3^2)^...
                %     (FUN_get_prop_NASA9(T3,'gamma')/(FUN_get_prop_NASA9(T3,'gamma')-1));
                ramp_prop(3,:) = [M3,theta3,beta3,P3,T3,P03];

                % check for pressure jump anlong with maxramp angle
                gamma = FUN_get_prop_NASA9(T3,'gamma');
                if theta(4) < FUN_oblique_shock(gamma, M3, 0, 'theta_max') && P3/P2 <= 1.7
                    % fourth shock
                    theta4 = theta(4);
                    [result] = FUN_oblique_shock_calculator_NASA9(theta4, M3, P3, T3, 'beta');
                    beta4 = result(1);
                    ratio1 = result(2:end-1);
                    M4 = result(end);
                    P4 = ratio1(1)*P3;
                    T4 = ratio1(2)*T3;
                    P04 = ratio1(3)*P03;
                    % P04 = P4 * (1 + 0.5*(FUN_get_prop_NASA9(T4,'gamma')-1)*M4^2)^...
                    %     (FUN_get_prop_NASA9(T4,'gamma')/(FUN_get_prop_NASA9(T4,'gamma')-1));
                    ramp_prop(4,:) = [M4,theta4,beta4,P4,T4,P04];

                    if P4/P3 < 1.8
                        % fifth shock
                        % have to straighten the flow
                        theta5 = theta1 + theta2 - theta3 - theta4;
                        gamma = FUN_get_prop_NASA9(T4,'gamma');
                        if theta5>0.1 && theta5<FUN_oblique_shock(gamma, M4, 0, 'theta_max')
                            [result] = FUN_oblique_shock_calculator_NASA9(theta5, M4, P4, T4, 'beta');
                            beta5 = result(1);
                            ratio1 = result(2:end-1);
                            M5 = result(end);
                            P5 = ratio1(1)*P4;
                            T5 = ratio1(2)*T4;
                            P05 = ratio1(3)*P04;
                            % P05 = P5 * (1 + 0.5*(FUN_get_prop_NASA9(T5,'gamma')-1)*M5^2)^...
                            %     (FUN_get_prop_NASA9(T5,'gamma')/(FUN_get_prop_NASA9(T5,'gamma')-1));
                            ramp_prop(5,:) = [M5,theta5,beta5,P5,T5,P05];

                            % check for pressure jump and temperature ratio
                            if P5/P4 < 1.8
                                [L,D,ramp_param, ramp_coord,coord_ramp1,Tth,Mth,wdt_in,wdt_cwl,lng_in, SI] = ...
                                    FUN_generate_ramp(alpha, M_oo, P_oo, T_oo, ramp_prop, h_th, m_dot, post,...
                                    config_param);

                                % Extracting PR, Tth, Mth, length, L/D, SI
                                % will be used ony for postprocess
                                postprocess = [ramp_param(5,4)/P_oo, Tth, Mth,...
                                    lng_in, L/D, SI];

                                if post=='y'  % this section used only during post processing
                                    fprintf('\nPR\t=\t%f',ramp_param(5,4)/P_oo)
                                    % fprintf('\nTR\t=\t%f',ramp_param(5,5)/T_oo)
                                    fprintf('\nTR\t=\t%f',Tth/T_oo)
                                    fprintf('\nTPR\t=\t%f',ramp_param(5,6)/P0_oo)
                                    fprintf('\nMth\t=\t%f',Mth)  %ramp_param(5,1))
                                    fprintf('\nTth (K)\t=\t%f',Tth)  %ramp_param(5,5))
                                    fprintf('\nDrag (N)\t=\t%f',D)
                                    fprintf('\nL/D\t=\t%f',L/D)
                                    fprintf('\nLength (m)\t=\t%f',lng_in)
                                    fprintf('\nWidth intake (m)\t=\t%f',wdt_in)
                                    fprintf('\nWidth cowl (m)\t=\t%f',wdt_cwl)
                                    fprintf('\nStartabililty index\t=\t%f\n\n',SI)
                                end

                                % constraint on length and throat temperature
                                % if (lng_in <= 4 && ramp_param(5,5)/T_oo >= 0.9*TR_th) ...
                                %         || post=='y'
                                if (Tth >= T_t && lng_in <= 4) ...
                                        || post=='y'
                                % if (lng_in <= 4) ...
                                %         || post=='y'
                                    obj_fn = [D, ramp_param(5,6)/P0_oo];  %store objective function
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

catch ME

    % catch any error in the internal functions (like loops exceeding max iterations) 
    % and design will be labelled as unfeasible
    % prevent crash of entire code
    % prevent loops in internal functions from running till an eternity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %while testing or adding new functions run unit tests
    %error will not be visible becasue of this
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!

    obj_fn = penalty;
    ramp_coord = 0;
    coord_ramp1 = 0;
    postprocess = [0, 0, 0, 0, 0, 0];
    return;
end

end