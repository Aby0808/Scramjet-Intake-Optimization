function [D] = FUN_compute_viscous_drag_force(gamma, M, P, T, Tw, l, S)

% this function calculates the skin friction force on the ramps
mu_r = 1.79*10^-5;
T_star = T*(1 + 0.032*M^2 + 0.58*(Tw/T) - 1);
mu_star = mu_r * ((T_star/293.15)^1.5 * (403.15/(T_star+110)));
rho_star = P/(287*T_star);

ue = M * sqrt(gamma*287*T);

Re_c = rho_star*ue*l/mu_star;

Cf = 1.328/sqrt(Re_c);

D = 0.5*Cf*rho_star*ue^2 *S; 
end