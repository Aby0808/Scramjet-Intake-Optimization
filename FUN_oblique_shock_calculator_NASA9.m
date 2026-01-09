function [value] = FUN_oblique_shock_calculator_NASA9(angle, M1, p1, T1, para)

% this function calculates the oblique shock properties from first principles

%% initializing properties

% Gas constant for air
Ru   = 8.31446261815324;       % J/mol-K
Mair = 28.966e-3;          % kg/mol
R    = Ru / Mair;

% === Upstream properties ===
gamma1  = FUN_get_prop_NASA9(T1, 'gamma');
a1      = sqrt(gamma1 * R * T1);
u1      = M1 * a1;

res = 100;
max_iter = 2000;
iter = 0;
ratio = [0 0 0 0 0 0 0];

%% solver

switch para
    case 'beta' % compute shock angle (theta) from deflection angle (beta)
         beta = FUN_oblique_shock(gamma1, M1, angle, 'beta');
         while res>10^-3 % find beta with this Newton_Raphson type solver
             u1n = u1*sind(beta);
             M1n = u1n/a1;
             u1t = u1*cosd(beta);
             [ratio, M2n] = FUN_normal_shock_calculator_NASA9(M1n, p1, T1);
             u2n = ratio(end)*u1n;
             u2t = u2n/tand(beta-angle);

             res = sqrt(((u2t-u1t)/u1t)^2); % beta such that tangential velocity matches upstream and downstream
             beta = beta+(u2t-u1t)/u2t;
             iter = iter + 1;
             if iter>max_iter
                 % fprintf('%f %f %f %f',angle, M1,p1,T1);
                 error('number of iterations went over the max iterations 2000');
                 % warning(['number of iterations went over the max iterations 2000. ' ...
                 %     'Reverting to oblique shock relation with constant gamma']);
                 % 
                 % % revert to oblique shock properties with constant gamma
                 % % to avoid code crashing at very low deflection angle
                 % gamma = FUN_get_prop_NASA9(T1,'gamma');
                 % beta = FUN_oblique_shock(gamma,M1,angle,'beta');
                 % M2 = FUN_oblique_shock(gamma,M1,angle,'M2');
                 % PR = FUN_oblique_shock_prop_ratio(gamma,M1,angle,'PR');
                 % TR = FUN_oblique_shock_prop_ratio(gamma,M1,angle,'TR');
                 % TPR = FUN_oblique_shock_prop_ratio(gamma,M1,angle,'TPR');
                 % DR = PR/TR;
                 % HR = (FUN_get_prop_NASA9(T1*TR,'cp')*T1*TR)/(FUN_get_prop_NASA9(T1,'cp')*T1);
                 % UR = (M2/sqrt(FUN_get_prop_NASA9(T1*TR,'gamma')*T1*TR))/...
                 %     (M1/sqrt(FUN_get_prop_NASA9(T1,'gamma')*T1));
                 % value = [beta, PR, TR, TPR, 1, DR, HR, UR, M2];
                 % return
             end
         end
         % finding post shock property ratios
         M2 = M2n/sind(beta-angle);
         % ratio(end) = ratio(end)/tand(beta-angle);
         value = [beta, ratio, M2];

    case 'theta' % Compute deflection angle (theta) from shock angle (beta)
        u1n = u1*sind(angle);
        M1n = u1n/a1;
        u1t = u1*cosd(angle);
        [ratio, M2n] = FUN_normal_shock_calculator_NASA9(M1n, p1, T1);
        u2n = ratio(end)*u1n;
        theta = angle - atand(u2n/u1t);
        M2 = M2n/sind(angle-theta);
        ratio(end) = ratio(end)/tand(angle-theta);
        value = [theta, ratio, M2];

    otherwise
        error('Invalid parameter requested: %s', parameter);
end

end