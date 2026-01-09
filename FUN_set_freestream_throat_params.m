function [param] = FUN_set_freestream_throat_params()

% this function sets the freestream and required throat parameters which
% will remain fixed for all calculations (for optimization and post processing)

% this function should be called in the beginning of main and postprocess program
% set params here

alpha = 4;
M_oo = 6.5;
P_oo = 1172;    %Pa
T_oo = 226.65;  %K
m_dot = 18.7;   %kg/s
h_th = 0.054;   %m
T_th = 1150;    %K

param = [alpha, M_oo, P_oo, T_oo, m_dot, h_th, T_th];

end