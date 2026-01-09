function [coord_intersect, coord_shock, T, M, mdot] = FUN_blunt_body_shock(R, Beta, Moo, Poo, Too, xcw, ycw, def, x, y, para)

% this function calculates the shape of blunt shock 
% it calculates the height of intake(x1,y1) or (x3,y3) intersection of ramp 2 and ramp 3 top 
% based on whether its a cowl or LE shock
% it also calculates the decrease in Mach number and increase in Temperature due to entropy layer

delta = @(R,M) R*(0.386*exp(4.67/M^2)); % shock stand-off distance
Rc = @(R,M) R*(1.386*exp(1.8/(M-1)^0.75)); % detached shock radius

switch para

    case 'cw'
        % case for cowl shock

        RotM = [cosd(def) -sind(def); sind(def) cosd(def)];  % rotation matrix to make cowl ramp horizontal
        RotMinv = transpose(RotM);  % rotation matrix to make intersection points back to roiginal coord system
        rot_coord = RotM*[x+xcw;y+ycw];
        xnew = -(R + delta(R,Moo) - Rc(R,Moo)*cotd(Beta)^2 * ((1 + (rot_coord(2)*tand(Beta)/Rc(R,Moo))^2)^0.5 - 1));
        coord_intersect = RotMinv*[xnew-xcw;rot_coord(2)-ycw];

        % find the bow shock curve
        ys = 0:0.001:rot_coord(2);
        xs = -(R + delta(R,Moo) - Rc(R,Moo)*cotd(Beta)^2 * ((1 + (ys*tand(Beta)/Rc(R,Moo)).^2).^0.5 - 1));

        coord_shock = RotMinv*[xs-xcw;ys-ycw];

    case 'le'
        % case for leading edge shock

        RotM = [cosd(-def) -sind(-def); sind(-def) cosd(-def)];  % rotation matrix to make cowl ramp horizontal
        RotMinv = transpose(RotM);  % rotation matrix to make intersection points back to roiginal coord system
        rot_coord = RotM*[xcw;ycw-y];
        xnew = -(R + delta(R,Moo) - Rc(R,Moo)*cotd(Beta)^2 * ((1 + (rot_coord(2)*tand(Beta)/Rc(R,Moo))^2)^0.5 - 1));
        coord_intersect = RotMinv*[-xnew;rot_coord(2)+y];

        % find the bow shock curve
        ys = 0:0.001:rot_coord(2);
        xs = -(R + delta(R,Moo) - Rc(R,Moo)*cotd(Beta)^2 * ((1 + (ys*tand(Beta)/Rc(R,Moo)).^2).^0.5 - 1));

        coord_shock = RotMinv*[xs;ys+y];

    otherwise
end

[T, M, mdot] = FUN_blunt_shock_entropy_layer(R, Beta, Moo, Poo, Too);

end