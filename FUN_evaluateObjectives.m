function obj = FUN_evaluateObjectives(param)

theta = param(1:end-1);
flag = param(end);
    % Constants
    gamma = 1.4;
    R = 287;
    M_oo = 6.5;
    P_oo = 1171;
    T_oo = 279;
    M_th = 2.1;
    m_dot = 18.7;
    PR_th = 100;

    % Call the objective function
    [obj_fn, coord, ~] = FUN_objective_function(gamma, R, M_oo, P_oo, T_oo, M_th, m_dot, PR_th, theta, flag, 'n');

    % if any(obj_fn==10^-6)
    %     ln_1 = 100;
    % else
    %     ln_1 = sqrt((coord(1,1) - coord(2,1))^2);
    % end

    % Extract objectives from FUN_objective_function
    drag = obj_fn(1);
    l2d = obj_fn(2);   % Lift-to-drag ratio
    pressureRatio = obj_fn(3);
    intakeLength = obj_fn(5); % Extracted from FUN_generate_ramp
    intakeWidth = obj_fn(6);  % Extracted from FUN_generate_ramp
    pressureRecovery = obj_fn(4); % Extracted from FUN_generate_ramp

    % Objective function (all to be minimized)
    obj = [drag, -l2d, -pressureRatio, intakeLength, intakeWidth, -pressureRecovery];
end
