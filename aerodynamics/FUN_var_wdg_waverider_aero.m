function [L,D,coord,Area,int_wdt] = FUN_var_wdg_waverider_aero(alpha, P_oo, ramp_param, ramp_coord, cwl_wdt)
%% Veriable Wedge Waverider Forebody Design
% based on the constant wedge waverider theory
%Starkey, Ryan P., and Mark J. Lewis. "Simple analytical model for parametric 
% studies of hypersonic waveriders." Journal of Spacecraft and Rockets 36.4 (1999): 516-523.

% Input Parametersin_wdt

theta = ramp_param(1,2);
length = sqrt((ramp_coord(1,1)-ramp_coord(2,1))^2); %+0.1    % Forebody length (m)

%% finding teh leading edge curve of the forebody
%leading edge is modelled using equation ax + b = y^4
%x1,y1 is the 

% X = [0 1; length 1];
% Y = [cwl_wdt^4; in_wdt^4];
% Coeff = X\Y;                                  % coefficients of the biquadratic eqn
% Coeff(2,1)=(cwl_wdt/2)^4;
% Coeff(1,1)=((in_wdt/2)^4 - Coeff(2,1))/length;
x1 = -0.49535*cwl_wdt/(2*tand(asind(1/ramp_param(1,1)))) + ramp_coord(2,1);
y1 = 1.49535*cwl_wdt/2;
x2 = ramp_coord(1,1); %-length;
y2 = cwl_wdt/2;
A = [x1 1; x2 1];
B = [y1^4 ; y2^4];
Coeff=inv(A)*B;

%% bottom compression surface
Area = 0;
dx=((Coeff(2)/Coeff(1)) + x1)/40; %0.05;
% Discretization
coord = zeros(1000,3);
c=1;
x=-(Coeff(2)/Coeff(1)):dx:x1;
for i=1:size(x,2)
    if i>1
        ymax = (Coeff(1)*x(i) + Coeff(2))^(1/4);
        yl = 0:ymax/(i*2):ymax;
    else
        ymax = 0;
        yl = 0;
    end

    if ~isreal(ymax)
        continue
    end
    
    if ~isreal(yl)
        disp('')
    end
    for j=1:size(yl,2)
        zl = (x(i)-ramp_coord(2,1))*tand(theta) - ramp_coord(2,2); %(ramp_coord(2,2)*(1-tand(theta)));
        coord(c,1) = x(i);
        coord(c,2) = yl(j);
        coord(c,3) = -zl;
        Area = Area + 2*ymax*dx/(i*10);
        c=c+1;
    end
end

if isempty(x)
    Area = (Area + (cwl_wdt*length))/cosd(theta);
    int_wdt = 0;
else
    Area = (Area + (0.5*(cwl_wdt+ymax)*length))/cosd(theta);
    int_wdt = 2*ymax;
end

force = (ramp_param(1,4) - P_oo)*Area;

force_visc = FUN_compute_viscous_drag_force(ramp_param(1,1),ramp_param(1,4),ramp_param(1,5),length,Area);

L = force*cosd(theta+alpha) - force_visc*sind(theta-alpha);
D = force*sind(theta+alpha) + force_visc*cosd(theta+alpha);

% mirror the y coordinate
coord=[coord; [coord(:,1)';-coord(:,2)';coord(:,3)']'];

% move the forbody to align it with the entire intake
% coord=[(coord(:,1)+ramp_coord(1,1))';(coord(:,3)+ramp_coord(1,2))';coord(:,2)']';
coord=[(coord(:,1))';(coord(:,3))';coord(:,2)']';

%% top surface

end