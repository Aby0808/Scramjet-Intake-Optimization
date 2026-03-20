function [M, P, T, P0] = FUN_compute_generalized_1d_flow(M1,P1,T1,l,D1,D2,A1,A2)

% this function computes the final state of flow properties due to area
% change, skin friction and wall heat flux

n=200;
dA = (A2-A1)/n;
dD = (D2-D1)/n;
dx = l/n;

gamma = FUN_get_prop_NASA9(T1,'gamma');
u = M1*sqrt(gamma*287*T1);
[P01, ~] = FUN_get_stagnation_properties(T1,P1,u);

M=M1; P=P1; T=T1; P0 = P01; A=A1; Dh = D1;

rho = P/(287*T);

for i = 1:n

    Tauw = FUN_compute_viscous_drag_force(M,P,T,Dh,1);
    f = 4*Tauw/(0.5*rho*u^2);

    dM2 = (-2*((1+0.5*(gamma-1)*M^2)/(1-M^2))*(dA/A) + ...
        (gamma*M^2 *(1+0.5*(gamma-1)*M^2)/(1-M^2))*(f*dx/Dh))*M^2;
    dP = (((gamma*M^2)/(1-M^2))*(dA/A) - ...
        (gamma*M^2 *(1+0.5*(gamma-1)*M^2)/(2*(1-M^2)))*(f*dx/Dh))*P;
    dT = ((((gamma-1)*M^2)/(1-M^2))*(dA/A) - ...
        (gamma*((gamma-1)*M^4)/(2*(1-M^2)))*(f*dx/Dh))*T;
    dP0 = -(0.5*gamma*M^2)*(f*dx/Dh)*P0;
    
    du = ((-1/(1-M^2))*(dA/A) + ((gamma*M^2)/(2*(1-M^2)))*(f*dx/Dh))*u;

    M = sqrt(M^2 + dM2);
    P = P + dP;
    T = T+dT;
    P0 = P0 + dP0;
    u = u + du;

    gamma = FUN_get_prop_NASA9(T,'gamma');
    Dh = Dh + dD;
    A = A + dA;

end

end