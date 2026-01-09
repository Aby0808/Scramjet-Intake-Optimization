function [value] = FUN_oblique_shock_calculator_avg_gamma(angle, M1, p1, T1, para)

% this function calculates the oblique shock properties from first principles

%% initializing properties

% Gas constant for air
Ru   = 8.31446261815324;       % J/mol-K
Mair = 28.966e-3;          % kg/mol
R    = Ru / Mair;

% === Upstream properties ===
gamma1  = FUN_get_prop_NASA9(T1, 'gamma');
cp1     = FUN_get_prop_NASA9(T1,'cp');
a1      = sqrt(gamma1 * R * T1);
u1      = M1 * a1;
h1      = cp1*T1;

%% solver

switch para
    case 'beta' % compute shock angle (theta) from deflection angle (beta)
        temp    = FUN_oblique_shock_calculator_NASA9(angle,M1,p1,T1,'beta');
        T2      = T1*temp(3);
        gamma2  = FUN_get_prop_NASA9(T2, 'gamma');
        % cp2     = FUN_get_prop_NASA9(T1,'cp');
        % a2      = sqrt(gamma2 * R * T2);
        gamma_avg = 1.4;%(gamma1+gamma2)/2;
        beta    = FUN_oblique_shock(gamma_avg, M1, angle, 'beta');
        % M2      = FUN_oblique_shock(gamma_avg, M1, angle, 'M2');
        % u2      = M2 * a2;
        % h2      = cp2*T2;
        % p2p1    = FUN_oblique_shock_prop_ratio(gamma_avg,M1,angle,'PR');
        % T2T1    = FUN_oblique_shock_prop_ratio(gamma_avg,M1,angle,'TR');
        ratio   = temp(2:end-1);
        value   = [beta, ratio, temp(end)];

    case 'theta' % Compute deflection angle (theta) from shock angle (beta)
        temp    = FUN_oblique_shock_calculator_NASA9(angle,M1,p1,T1,'theta');
        T2      = T1*temp(3);
        gamma2  = FUN_get_prop_NASA9(T2, 'gamma');
        % cp2     = FUN_get_prop_NASA9(T1,'cp');
        % a2      = sqrt(gamma2 * R * T2);
        gamma_avg = 1.4; %(gamma1+gamma2)/2;
        theta    = FUN_oblique_shock(gamma_avg, M1, angle, 'theta');
        % M2      = FUN_oblique_shock(gamma_avg, M1, theta, 'beta2M2');
        % u2      = M2 * a2;
        % h2      = cp2*T2;
        % p2p1    = FUN_oblique_shock_prop_ratio(gamma_avg,M1,theta,'PR');
        % T2T1    = FUN_oblique_shock_prop_ratio(gamma_avg,M1,theta,'TR');
        ratio   = temp(2:end-1);
        value   = [theta, ratio, temp(end)];

    otherwise
        error('Invalid parameter requested: %s', parameter);
end

end