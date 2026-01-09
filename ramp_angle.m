clear
clc
close all

%% input parameters
gamma = 1.4;

M1 = 6.5;
P1 = 1197;
T1 = 226.5;

PR = 80;
Tmax = 1000;
mdot = 5;

%% calculating the shock angles

for beta1 = asind(1/M1):0.1:FUN_oblique_shock(gamma, M1, 0, 'beta_max')
    
end