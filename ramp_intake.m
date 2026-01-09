clear
clc
close all

%% 

% this code finds the ramp angles for intake

%% initialization

Moo = 6.5;
Poo = 1017;
Too = 279;
mdot = 15;
gamma = 1.4;
alpha = 4;
Mt = 2;
Tt = 900;

data = zeros(10,100);

%% interating to find ramp angles

theta = 0.01:0.01:15;

for i=1:size(theta,2)
    M2 = FUN_oblique_shock(gamma,Moo,theta(i),'M2');
    for j = 1:size(theta,2)
        
    end
end