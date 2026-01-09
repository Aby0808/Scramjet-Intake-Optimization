function [value] = FUN_conical_shock_prop_ratio(gamma, M, theta, parameter)
% FUN_conical_shock_prop_ratio: Computes pressure, temperature, and total pressure ratios across a conical shock.

% Compute beta (shock angle) using the Taylor-Maccoll solution
beta = FUN_conical_shock(gamma, M, theta, 'beta');

% Compute normal Mach number (Mn1) using beta
Mn1 = M * sind(beta);

Mc = FUN_conical_shock(gamma,M,theta,'theta2M2'); %Mach at cone surface
M_ob = FUN_oblique_shock(gamma,M,beta,'beta2M2'); %Mach post shock

switch parameter
    case 'PR'
        PR_ob = (2*gamma*Mn1^2 - (gamma-1))/(gamma+1); %PR accross shock
        value = PR_ob*((1+0.5*(gamma-1)*M_ob^2)/(1+0.5*(gamma-1)*Mc^2)) ...
            ^(gamma/(gamma-1));

    case 'TR'
        % Static temperature ratio across the normal shock
        TR_ob = (2*gamma*Mn1^2 - (gamma-1)) * ((gamma-1)*Mn1^2 + 2) ...
            / ((gamma+1)^2 * Mn1^2);
        value = TR_ob*((1+0.5*(gamma-1)*M_ob^2)/(1+0.5*(gamma-1)*Mc^2));

    case 'TPR'
        % Total pressure ratio across the shock
        PR_ob = (2*gamma*Mn1^2 - (gamma-1))/(gamma+1); %PR accross shock
        value = PR_ob*((1+0.5*(gamma-1)*M_ob^2)/(1+0.5*(gamma-1)*M^2)) ...
            ^(gamma/(gamma-1));

    otherwise
        error('Wrong parameter specified: %s. Use ''PR'', ''TR'', or ''TPR''.', parameter);
end

end
