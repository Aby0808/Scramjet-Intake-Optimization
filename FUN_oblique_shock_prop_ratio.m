function [value] = FUN_oblique_shock_prop_ratio(gamma, M, theta, parameter)
% oblique_shock_ratios: Computes the pressure and temperature ratios across an oblique shock.

% Compute beta (shock angle) using the existing function
beta = FUN_oblique_shock(gamma, M, theta, 'beta');

% Compute normal Mach number (Mn1)
Mn1 = M * sind(beta);

switch parameter
    case 'PR'
        % Pressure ratio across normal shock
        value = (2*gamma*Mn1^2 - (gamma-1))/(gamma+1);

    case 'TR'
        % Temperature ratio across normal shock
        value = (2*gamma*Mn1^2 - (gamma-1))*((gamma-1)*Mn1^2 + 2) ...
            / ((gamma+1)^2 *Mn1^2);
    case 'TPR'
         % Total pressure ratio across the shock
         value = FUN_oblique_shock_prop_ratio(gamma, M, theta, 'PR') * ...
             ((1 + 0.5*(gamma-1)*FUN_oblique_shock(gamma, M, theta, 'theta2M2')^2) ...
             /(1 + 0.5*(gamma-1)*M^2))^(gamma/(gamma-1));
    otherwise
        error('\n\n wrong parameter: %s specified',parameter)
end
