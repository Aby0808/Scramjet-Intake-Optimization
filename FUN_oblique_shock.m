function [value] = FUN_oblique_shock(gamma, M, angle, parameter)
% Robust oblique shock solver using the bisection method for weak shock solutions

switch parameter
    case 'beta'
        % Convert deflection angle to radians
        theta = deg2rad(angle);

        % Check if deflection angle is within physical limits
        theta_max = FUN_oblique_shock(gamma, M, 0, 'theta_max');
        if angle > theta_max || angle < 0
            error('Deflection angle is out of bounds for the given Mach number.');
        end

        % Define bounds for weak shock solution
        beta_min = asin(1 / M); % Mach angle (lower bound)
        beta_max = pi / 2 - 1e-6; % Slightly less than 90 degrees (upper bound)

        % Function for f(beta)
        f_beta = @(beta) 2 * cot(beta) * ((M^2 * sin(beta)^2 - 1) / ...
            (M^2 * (gamma + cos(2 * beta)) + 2)) - tan(theta);

        % Use fsolve with an initial guess near the Mach angle
        options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10);
        beta_guess = beta_min + 0.1; % Slightly above Mach angle

        % Solve for beta using fsolve
        beta_solution = fsolve(f_beta, beta_guess, options);

        % Return beta in degrees
        value = rad2deg(beta_solution);

    % case 'beta' % Compute shock angle (beta) from deflection angle (theta)
    %     theta = deg2rad(angle); % Convert deflection angle to radians
    % 
    %     % Check if deflection angle is within physical limits
    %     theta_max = FUN_oblique_shock(gamma, M, 0, 'theta_max');
    %     if angle > theta_max || angle < 0
    %         error('Deflection angle is out of bounds for the given Mach number.');
    %     end
    % 
    %     % Define bounds for beta (weak shock solution)
    %     beta_min = asin(1 / M); % Mach angle (lower bound)
    %     beta_max = pi / 2 - 1e-6; % Slightly less than 90 degrees (upper bound)
    % 
    %     % Function for f(beta)
    %     f_beta = @(beta) 2 * cot(beta) * ((M^2 * sin(beta)^2 - 1) / ...
    %         (M^2 * (gamma + cos(2 * beta)) + 2)) - tan(theta);
    % 
    %     % Bisection method for solving f(beta) = 0
    %     for i = 1:500 % Increased max iterations
    %         beta_mid = (beta_min + beta_max) / 2; % Midpoint of the range
    %         f_mid = f_beta(beta_mid);
    % 
    %         % Convergence check
    %         if abs(f_mid) < 1e-10 % Tightened tolerance
    %             value = rad2deg(beta_mid); % Convert to degrees
    %             return;
    %         end
    % 
    %         % Update bounds
    %         if sign(f_mid) == sign(f_beta(beta_min))
    %             beta_min = beta_mid;
    %         else
    %             beta_max = beta_mid;
    %         end
    %     end
    % 
    %     % If no convergence
    %     error('Bisection method did not converge after 500 iterations.');


    case 'theta' % Compute deflection angle (theta) from shock angle (beta)
        beta = deg2rad(angle); % Convert shock angle to radians
        value = atand(2 * cot(beta) * ((M^2 * sin(beta)^2 - 1) / ...
                   (M^2 * (gamma + cos(2 * beta)) + 2)));

    case 'M2' % Compute post-shock Mach number (M2)
        beta = FUN_oblique_shock(gamma, M, angle, 'beta'); % Get shock angle
        Mn1 = M * sind(beta); % Normal Mach number before the shock
        Mn2 = sqrt((1 + 0.5 * (gamma - 1) * Mn1^2) / ...
                   (gamma * Mn1^2 - (gamma - 1) * 0.5)); % Normal Mach number after the shock
        value = Mn2 / sind(beta - deg2rad(angle)); % Convert back to oblique Mach number

    case 'theta2M2' % Compute post-shock Mach number (M2) from deflection angle (theta)
        beta = FUN_oblique_shock(gamma, M, angle, 'beta'); % Get shock angle
        Mn1 = M * sind(beta);
        Mn2 = sqrt((1 + 0.5 * (gamma - 1) * Mn1^2) / (gamma * Mn1^2 - (gamma - 1) * 0.5));
        value = Mn2 / sind(beta - angle);

    case 'beta2M2' % Compute post-shock Mach number (M2) from wave angle (beta)
        beta = angle;
        theta = FUN_oblique_shock(gamma, M, angle, 'theta'); % Get deflection angle
        Mn1 = M * sind(beta);
        Mn2 = sqrt((1 + 0.5 * (gamma - 1) * Mn1^2) / (gamma * Mn1^2 - (gamma - 1) * 0.5));
        value = Mn2 / sind(beta - theta);

    case 'beta_max' % Maximum shock angle (beta_max)
        beta_min = asin(1 / M); % Mach angle
        beta_max = pi / 2 - 1e-6; % Slightly less than 90 degrees
        beta_guess = linspace(beta_min, beta_max, 1000);
        theta_values = arrayfun(@(beta) atand(2 * cot(beta) * ...
            ((M^2 * sin(beta)^2 - 1) / (M^2 * (gamma + cos(2 * beta)) + 2))), beta_guess);
        [~, max_idx] = max(theta_values); % Find maximum deflection angle
        value = rad2deg(beta_guess(max_idx));

    case 'theta_max' % Maximum deflection angle (theta_max)
        beta_min = asin(1 / M); % Mach angle
        beta_max = pi / 2 - 1e-6; % Slightly less than 90 degrees
        beta_guess = linspace(beta_min, beta_max, 1000);
        theta_values = arrayfun(@(beta) atand(2 * cot(beta) * ...
            ((M^2 * sin(beta)^2 - 1) / (M^2 * (gamma + cos(2 * beta)) + 2))), beta_guess);
        value = max(theta_values);

    otherwise
        error('Invalid parameter requested: %s', parameter);
end
end
