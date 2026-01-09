function val = FUN_get_prop_NASA9(T,para)

% this function calculates Cp from NASA9 polynomials
% air_thermo_nasa7  Return cp [J/kg-K] or h [J/kg] for AIR (NASA-7).
% Inputs:  T [K] (scalar or vector), what = 'cp' or 'h'
% Range:   200–6000 K (piecewise at 1000 K). Outside that, it extrapolates.
% Data:    Gordon (1982) mix: Air (N2 78.084, O2 20.9476, Ar 0.9365, CO2 0.0319)
% M = 28.9651159 g/mol

Ru   = 8.31446261815324;     % J/mol-K
M    = 28.9651159e-3;        % kg/mol
Rsp  = Ru/M;                 % J/kg-K
Tmid = 1000.0;

% NASA 9 coefficient format (NASA Glenn)
aL = [ 1.009950160e+04, -1.968275610e+02,  5.009155110e+00, ...
      -5.761013730e-03,  1.066859930e-05, -7.940297970e-09, ...
      2.185231910e-12, -1.76796731e+2, -3.921504225 ];
aH = [ 2.415214430e+05, -1.257874600e+03,  5.144558670e+00, ...
      -2.138541790e-04,  7.065227840e-08, -1.071483490e-11, ...
      6.577800150e-16, 6.46226319e+3, -8.147411905];

% Chemkin 7 coefficient format (Gordon McBride)
% aL = [3.57e+00, -6.79e-04, 1.55E-06, -3.30e-12, -4.66e-13, ...
%     -1.06e+03, 3.72e+00];
% aH = [3.09e+00, 1.25e-03, -4.24e-07, 6.75e-11, -3.97e-15, ...
%     -9.95e+02, 5.96e+00];

% Coefficients taken from Fluent
% aL = [2898903, -56496.26, 1437.799, -1.653609, 0.003062254, ...
%     -2.279138e-6, 6.272365e-10];
% aH = [6.932494e+7, -361053.2, 1476.665, -0.06138349, 2.027963e-5, ...
%     -3.075525e-9, 1.888054e-13];

Tcol = T(:);
useL = (Tcol <= Tmid);
useH = ~useL;

% NASA 9 coefficient format (NASA Glenn)
cp_over_R = @(a,TT) a(1).*TT.^(-2) + a(2).*TT.^(-1) + a(3) + a(4).*TT ...
                  + a(5).*TT.^2 + a(6).*TT.^3 + a(7).*TT.^4;

h_over_RT = @(a,TT) -a(1).*TT.^(-2) + a(2).*log(TT)./TT + a(3) ...
                  + a(4).*TT./2 + a(5).*TT.^2./3 + a(6).*TT.^3./4 ...
                  + a(7).*TT.^4./5 + a(8)./TT;

s_over_RT = @(a,TT) -0.5*a(1).*TT.^(-2) - a(2).*TT.^(-1) + a(3).*log(TT) ...
                  + a(4).*TT + a(5).*TT.^2./2 + a(6).*TT.^3./3 ...
                  + a(7).*TT.^4./4 + a(9);

% Chemkin 7 coefficient format (Gordon McBride)
% cp_over_R = @(a,T) a(1) + a(2)*T + a(3)*T^2 + a(4)*T^3 + a(5)*T^4;
% 
% h_over_RT = @(a,T) a(1) + a(2)*T/2 + a(3)*T^2 /3 + a(4)*T^3 /4 + ...
%     (a(5)*T^4)/5 + a(6)/T;

% Fluent polynomial
% cp_over_R = @(a,T) a(1)*T^6 + a(2)*T^5 + a(3)*T^4 + a(4)*T^3 + a(5)*T^2 ...
%      + a(6)*T + a(7);
% h_over_RT = @(a,T) a(1)*T^7 + a(2)*T^6 + a(3)*T^5 + a(4)*T^4 + a(5)*T^3 ...
%      + a(6)*T^2 + a(7)*T;


val = zeros(size(Tcol));
switch lower(para)
    case 'cp'
        if any(useL), val(useL) = Rsp .* cp_over_R(aL, Tcol(useL)); end
        if any(useH), val(useH) = Rsp .* cp_over_R(aH, Tcol(useH)); end
    case 'h'
        if any(useL), val(useL) = Rsp .* Tcol(useL) .* h_over_RT(aL, Tcol(useL)) + 302148; end  %absolute enthalpy
        if any(useH), val(useH) = Rsp .* Tcol(useH) .* h_over_RT(aH, Tcol(useH)) + 302148; end  %enthalpy of formation = 302148J/Kg
    case 's'
        if any(useL), val(useL) = Rsp .* s_over_RT(aL, Tcol(useL)); end
        if any(useH), val(useH) = Rsp .* s_over_RT(aH, Tcol(useH)); end
    case 'gamma'
        cp = FUN_get_prop_NASA9(T,'cp');
        val = cp/(cp - Rsp);
    case 'k'
        val = 0.05;  % assumed constant for now
    otherwise
        error("Argument must be 'cp', 'h', 'gamma' or 'h'.");
end

% restore input shape
if isrow(T), val = val.';
end

end