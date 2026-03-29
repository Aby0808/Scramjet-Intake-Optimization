function param = FUN_set_freestream_throat_params()

% this function sets the freestream and required throat parameters which
% will remain fixed for all calculations (for optimization and post processing)

% this function should be called in the beginning of main and postprocess program
% set params here

alpha = 2;
M_oo = 7.0;
P_oo = 1172;    %Pa
T_oo = 226.65;  %K
m_dot = 10;   %kg/s
h_th = 0.06;   %m
T_th = 1000;    %K

param = struct( ...
    'alpha', alpha, ...
    'M_oo', M_oo, ...
    'P_oo', P_oo, ...
    'T_oo', T_oo, ...
    'm_dot', m_dot, ...
    'h_th', h_th, ...
    'T_th', T_th, ...
    'vector', [alpha, M_oo, P_oo, T_oo, m_dot, h_th, T_th]);

end
