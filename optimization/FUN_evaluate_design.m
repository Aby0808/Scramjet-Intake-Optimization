function result = FUN_evaluate_design(alpha, M_oo, P_oo, T_oo, m_dot, h_th, T_t, theta, post, config_param)

% Full design-evaluation pipeline used by both optimization and post-processing.
% This keeps the current logic intact while separating design analysis from wrappers.
% The sequencing here follows the original objective-function logic:
% 1. check ramp-angle feasibility
% 2. march the flow through the shock system
% 3. generate the intake geometry
% 4. compute post-processing metrics
% 5. apply hard constraints and package the final result

%% initializing required parameters
% TR_th = T_t/T_oo;
gamma = FUN_get_prop_NASA9(T_oo,'gamma');

% Compute freestream properties
u_oo = M_oo*sqrt(FUN_get_prop_NASA9(T_oo,'gamma')*287*T_oo);
[P0_oo, ~] = FUN_get_stagnation_properties(T_oo,P_oo,u_oo);
% P0_oo = P_oo * (1 + 0.5*(gamma-1)*M_oo^2)^(gamma/(gamma-1));  % for some reason the
% stagnation properties are not coming properly

penalty = [10^5, -10^5];  % penalty to be returned for undesired flowfield

result = buildFailureResult(penalty);
ramp_prop = zeros(5,6);

try
    %% ramp generation and filtering

    % net deflection from first two ramp should be greater than deflection of third ramp
    if theta(1)+theta(2)-theta(3) <= 0
        result = buildFailureResult(penalty, ...
            "Net deflection after ramp 3 is non-positive.", "initial_ramp_feasibility");
        return
    end

    if theta(1)+theta(2)-theta(3)-theta(4) <= 0
        result = buildFailureResult(penalty, ...
            "Straightening ramp deflection becomes non-positive.", "initial_ramp_feasibility");
        return
    end

    if theta(1) >= FUN_oblique_shock(gamma, M_oo, 0, 'theta_max')
        result = buildFailureResult(penalty, ...
            "Ramp 1 exceeds the maximum oblique-shock turning angle.", "shock_1_filter");
        return
    end

    % first ramp
    theta1 = theta(1);   % preserve shock strength and align beta back to coord system
    resultShock = FUN_oblique_shock_calculator_NASA9(theta1+alpha, M_oo, P_oo, T_oo, 'beta');
    beta1 = resultShock(1) - alpha;
    ratio1 = resultShock(2:end-1);
    M1 = resultShock(end);
    P1 = ratio1(1)*P_oo;
    T1 = ratio1(2)*T_oo;
    P01 = ratio1(3)*P0_oo;
    ramp_prop(1,:) = [M1,theta1,beta1,P1,T1,P01];

    % second ramp
    % filtering based on max ramp angle
    gamma = FUN_get_prop_NASA9(T1,'gamma');
    if theta(2) >= FUN_oblique_shock(gamma, M1, 0, 'theta_max')
        result = buildFailureResult(penalty, ...
            "Ramp 2 exceeds the maximum oblique-shock turning angle.", "shock_2_filter");
        return
    end

    theta2 = theta(2);
    resultShock = FUN_oblique_shock_calculator_NASA9(theta2, M1, P1, T1, 'beta');
    beta2 = resultShock(1);
    ratio1 = resultShock(2:end-1);
    M2 = resultShock(end);
    P2 = ratio1(1)*P1;
    T2 = ratio1(2)*T1;
    P02 = ratio1(3)*P01;
    ramp_prop(2,:) = [M2,theta2,beta2,P2,T2,P02];

    % third shock
    gamma = FUN_get_prop_NASA9(T2,'gamma');
    if theta(3) >= FUN_oblique_shock(gamma, M2, 0, 'theta_max')
        result = buildFailureResult(penalty, ...
            "Ramp 3 exceeds the maximum oblique-shock turning angle.", "shock_3_filter");
        return
    end

    theta3 = theta(3);
    resultShock = FUN_oblique_shock_calculator_NASA9(theta3, M2, P2, T2, 'beta');
    beta3 = resultShock(1);
    ratio1 = resultShock(2:end-1);
    M3 = resultShock(end);
    P3 = ratio1(1)*P2;
    T3 = ratio1(2)*T2;
    P03 = ratio1(3)*P02;
    ramp_prop(3,:) = [M3,theta3,beta3,P3,T3,P03];

    % fourth shock
    % check for pressure jump along with max ramp angle
    gamma = FUN_get_prop_NASA9(T3,'gamma');
    if theta(4) >= FUN_oblique_shock(gamma, M3, 0, 'theta_max')
        result = buildFailureResult(penalty, ...
            "Ramp 4 exceeds the maximum oblique-shock turning angle.", "shock_4_filter");
        return
    end

    if P3/P2 > 1.7
        result = buildFailureResult(penalty, ...
            "Pressure jump across shock 3 exceeds the allowed limit.", "shock_4_filter");
        return
    end

    theta4 = theta(4);
    resultShock = FUN_oblique_shock_calculator_NASA9(theta4, M3, P3, T3, 'beta');
    beta4 = resultShock(1);
    ratio1 = resultShock(2:end-1);
    M4 = resultShock(end);
    P4 = ratio1(1)*P3;
    T4 = ratio1(2)*T3;
    P04 = ratio1(3)*P03;
    ramp_prop(4,:) = [M4,theta4,beta4,P4,T4,P04];

    if P4/P3 >= 1.8
        result = buildFailureResult(penalty, ...
            "Pressure jump across shock 4 exceeds the allowed limit.", "shock_4_filter");
        return
    end

    % fifth shock
    % have to straighten the flow
    theta5 = theta1 + theta2 - theta3 - theta4;
    gamma = FUN_get_prop_NASA9(T4,'gamma');
    if theta5 <= 0.1
        result = buildFailureResult(penalty, ...
            "Straightening ramp angle is too small.", "shock_5_filter");
        return
    end

    if theta5 >= FUN_oblique_shock(gamma, M4, 0, 'theta_max')
        result = buildFailureResult(penalty, ...
            "Straightening ramp exceeds the maximum oblique-shock turning angle.", "shock_5_filter");
        return
    end

    resultShock = FUN_oblique_shock_calculator_NASA9(theta5, M4, P4, T4, 'beta');
    beta5 = resultShock(1);
    ratio1 = resultShock(2:end-1);
    M5 = resultShock(end);
    P5 = ratio1(1)*P4;
    T5 = ratio1(2)*T4;
    P05 = ratio1(3)*P04;
    ramp_prop(5,:) = [M5,theta5,beta5,P5,T5,P05];

    % check for pressure jump and temperature ratio
    if P5/P4 >= 1.8
        result = buildFailureResult(penalty, ...
            "Pressure jump across shock 5 exceeds the allowed limit.", "shock_5_filter");
        return
    end

    %% geometry generation, viscous/blunt corrections, and aero forces
    [L,D,ramp_param, ramp_coord,coord_ramp1,Tth,Mth,wdt_in,wdt_cwl,lng_in, SI] = ...
        FUN_generate_ramp(alpha, M_oo, P_oo, T_oo, ramp_prop, h_th, m_dot, post, config_param);

    % Extracting PR, Tth, Mth, length, L/D, SI
    % will be used only for postprocess
    postprocess = [ramp_param(5,4)/P_oo, Tth, Mth, lng_in, L/D, SI];

    %% reporting block for post-processing runs
    if post=='y'
        FUN_print_design_summary(P_oo, T_oo, P0_oo, ramp_param, Tth, Mth, D, L, lng_in, wdt_in, wdt_cwl, SI)
    end

    %% hard constraints and final objective packaging
    % constraint on length and throat temperature
    % if (lng_in <= 4 && ramp_param(5,5)/T_oo >= 0.9*TR_th) || post=='y'
    if (Tth >= T_t && lng_in <= 4) || post=='y'
        result = struct( ...
            'objectiveFunction', [D, ramp_param(5,6)/P0_oo], ...
            'rampCoord', ramp_coord, ...
            'coordRamp1', coord_ramp1, ...
            'postprocess', postprocess, ...
            'isFeasible', true, ...
            'failureReason', "", ...
            'failureStage', "", ...
            'theta', theta, ...
            'rampProperties', ramp_param);
    else
        result = buildFailureResult(penalty, ...
            "Hard constraints on throat temperature or intake length were not satisfied.", ...
            "final_constraints");
    end

catch ME %#ok<NASGU>
    % catch any error in the internal functions (like loops exceeding max iterations)
    % and the design will be labelled as unfeasible
    % prevent crash of entire code
    % prevent loops in internal functions from running till an eternity
    %
    % while testing or adding new functions run unit tests
    % error will not be visible because of this
    result = buildFailureResult(penalty, ...
        "Unhandled internal exception during design evaluation.", "exception");
end

end

function result = buildFailureResult(penalty, failureReason, failureStage)
if nargin < 2
    failureReason = "Infeasible design.";
end

if nargin < 3
    failureStage = "";
end

result = struct( ...
    'objectiveFunction', penalty, ...
    'rampCoord', 0, ...
    'coordRamp1', 0, ...
    'postprocess', [0, 0, 0, 0, 0, 0], ...
    'isFeasible', false, ...
    'failureReason', failureReason, ...
    'failureStage', failureStage, ...
    'theta', [], ...
    'rampProperties', []);
end
