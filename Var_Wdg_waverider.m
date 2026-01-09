%close all

%% Constant Wedge Waverider Forebody Design
% based on the constant wedge waverider theory
%Starkey, Ryan P., and Mark J. Lewis. "Simple analytical model for parametric 
% studies of hypersonic waveriders." Journal of Spacecraft and Rockets 36.4 (1999): 516-523.

% Input Parameters
theta_deg = -15;  % Wedge angle in degrees
M = 6.5;        % freestream mach number
n = 0.5;        % Power-law exponent for planform
m = 0.5;        % power law exponent for variable wedge
length = 1.28;    % Vehicle length (m)
width_t = 2*0.631333;     % Vehicle width (m)
width_i = 2*0.294028;
delta = ramp_theta(1);      % Leading edge angle in degrees (for variable wedge) (value)
gamma = 1.4;
type = 'cc';    % convex cv / concave cc
width_w = width_t - width_i;
%% solving oblique shock relations ot find beta

beta = FUN_oblique_shock(gamma, M, -theta_deg, 'beta');

%%
theta = theta_deg;
A = width_w / (2 * length^n);  % Power-law constant for planform
B = A/tand(beta)^n;

% Discretization
coord = zeros(1000,3);
c=2;
x=0.0:0.02:length;
for i=1:size(x,2)
    switch type
        case 'cv'
            C = (A*x(i)^(n-m))/sqrt((tand(theta) - tand(delta))^(2*m));
        case 'cc'
            C = (A*x(i)^(n-m))/sqrt((tand(theta) + tand(-delta))^(2*m));
        otherwise
            error('wrong type selected: %s',type)
    end

    ymax = width_i/2 + A*x(i)^n;
    yl = 0:0.01:width_i/2;
    % flat portion of waverider
    for j=1:size(yl,2)
        zl = x(i)*tand(theta);
        coord(c,1) = x(i);
        coord(c,2) = yl(j);
        coord(c,3) = -zl;
        c=c+1;
    end

    yl = width_i/2:ymax/(i):ymax;
    for j=1:size(yl,2)
        zl = x(i)*tand(theta) - (yl(j)/C)^(1/m);
        coord(c,1) = x(i);
        coord(c,2) = yl(j);
        coord(c,3) = -zl;
        c=c+1;
    end
end

% Create Surface Plot Data
figure;
scatter3(coord(:,1),coord(:,2),coord(:,3))
hold on
scatter3(coord(:,1),-coord(:,2),coord(:,3))
% Adjust View and Labels
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
grid on;
axis equal;
% hold off;
