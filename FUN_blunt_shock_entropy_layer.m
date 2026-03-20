function [T,M,mdot] = FUN_blunt_shock_entropy_layer(R, Beta, Moo, Poo, Too)

% this function calculates the increase in temperature and decrease in mach
% number mass averaged

delta = @(R,M) R*(0.386*exp(4.67/M^2)); % shock stand-off distance
Rc = @(R,M) R*(1.386*exp(1.8/(M-1)^0.75)); % detached shock radius

y = 0:0.01:0.2;

x = -(R + delta(R,Moo) - Rc(R,Moo)*cotd(Beta)^2 * ((1 + (y.*tand(Beta)/Rc(R,Moo)).^2).^0.5 - 1));
M = 0;
T = 0;
mdot = 0;
for i=2:size(y,2)-1
    slope = (y(i)-y(i-1))/(x(i)-x(i-1));
    angle = atand(slope);
    Mn = Moo*cosd(angle);
    uoo = Mn*sqrt(FUN_get_prop_NASA9(Too,'gamma')*287*Too);
    rhooo = Poo/(287*Too);
    % [ratio, M2] = FUN_normal_shock_calculator_NASA9(Mn,Poo,Too);
    % T2 = ratio(2)*Too;
    % u2 = uoo * ratio(end);
    % rho2 = rhooo * ratio(end-2);

    % using constant gamma normal shock relation instead of variable gamma
    % using variable gamma significantly increased execution time without gaining significant accuracy
    u2 = uoo * (2 + 0.4*Mn^2)/(2.4*Mn^2);
    rho2 = rhooo * (2.4*Mn^2)/(2 + 0.4*Mn^2);
    T2 = Too * (1 + 2.8*(Mn^2  - 1)/2.4) * ((2 + 0.4*Mn^2)/(2.4*Mn^2));
    M2 = u2 / sqrt(1.4*287*T2);

    Ar = sqrt((y(i)-y(i-1))^2 + (x(i)-x(i-1))^2);
    mdot = mdot + rho2*u2*Ar;
    T = T + T2*rho2*u2*Ar;
    M = M + M2*rho2*u2*Ar;
    if angle < 1.5*Beta
        break;
    end
end

T = T/mdot;
M = M/mdot;

end