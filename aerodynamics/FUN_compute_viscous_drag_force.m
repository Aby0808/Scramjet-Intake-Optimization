function [D] = FUN_compute_viscous_drag_force(M, P, T, l, S)

% this function calculates the skin friction force on the ramps
gamma = FUN_get_prop_NASA9(T,'gamma');
Cp = FUN_get_prop_NASA9(T,'cp');
k = FUN_get_prop_NASA9(T,'k');
mu_r = 1.79*10^-5;

Pr = mu_r*Cp/k;
r = Pr^(1/3);  % recovery factor

u = M*sqrt(gamma*287*T);
hw = FUN_get_prop_NASA9(T,'h') + r*0.5*u^2;
Tw = FUN_get_T_from_h(hw);
% Tw = T*(1 + r*0.5*(gamma-1)*M^2);  % recovery wall temperature for adiabatic walls

T_star = T*(1 + 0.032*M^2 + 0.58*(Tw/T) - 1);
mu_star = mu_r * ((T_star/293.15)^1.5 * (403.15/(T_star+110))); %sutherland
rho_star = P/(287*T_star);

ue = M * sqrt(gamma*287*T);

Re_c = rho_star*ue*l/mu_star;

Cf = 1.328/sqrt(Re_c);

D = 0.5*Cf*rho_star*ue^2 *S; 

end