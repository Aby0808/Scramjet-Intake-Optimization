function [ratio, M2] = FUN_normal_shock_calculator_NASA9(M1, p1, T1)

% This function solves the normal shock properties using first principles

%% initializing properties

% Gas constant for air
Ru   = 8.31446261815324;       % J/mol-K
Mair = 28.966e-3;          % kg/mol
R    = Ru / Mair;

% === Upstream properties ===
gamma1  = FUN_get_prop_NASA9(T1, 'gamma');
a1      = sqrt(gamma1 * R * T1);
u1      = M1 * a1;
rho1    = p1 / (R * T1);
h1      = FUN_get_prop_NASA9(T1, 'h');
[P01,T01] = FUN_get_stagnation_properties(T1,p1,u1);

T2 = T1;
rho12 = (((gamma1+1)*M1^2)/(2 + (gamma1-1)*M1^2))^-1;
res = 100; 
rho_res = rho12;
max_iter = 2000;
iter = 0;
%% solver

while res>10^-3
    p2 = p1 + rho1*u1^2 *(1 - rho12);
    h2 = h1 + (u1^2 / 2)*(1 - rho12^2);

    T2 = FUN_get_T_from_h(h2);

    rho2 = p2/(T2*R);
    rho12 = rho1/rho2;
    res = sqrt(((rho12-rho_res)/rho12)^2);

    iter = iter + 1;
    if iter>max_iter
        error('number of iterations went over the max iterations 1000');
    end
    rho_res = rho12;
end

%% results

gamma2  = FUN_get_prop_NASA9(T2, 'gamma');
u2      = rho1*u1/rho2;
a2      = sqrt(gamma2 * R * T2);
M2      = u2 / a2;
[P02,T02] = FUN_get_stagnation_properties(T2,p2,u2);
ratio   = [p2/p1, T2/T1, P02/P01, T02/T01, rho2/rho1, h2/h1, u2/u1];

end