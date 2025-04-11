clear all; close all; clc;

%% Parameter Definitions (Scoring Function)

rWp = 2.25*10^(-1); %Revenue per unit payload weight
rVp = 3.00*10^2; %Revenue per unit payload volume
c_e = 8.33*10^(-8); %Cost per unit energy
c_c = 5.00*10^(-1); %Cost (COC) per flight
cWg = 6.74*10^1; %Cost per unit gross weight
c_f = 5.56*10^(-5); %Cost (FOC) per unit flight-time
t_l = 3.60*10^6; %Total flight time in aircraft life

%% Chosen Variables
L2D = 12; %Lift to drag ratio
EWF = 0.65; %Empty weight fraction
nb = 0.93; %efficiency
v = 20; %velocity

%% Establish Function

J = @(Wg,Ef) Wg*(nb*v*((rWp*Wg*(1-EWF)+rVp*Wg^(3/2)-c_e*Ef-c_c)/(Ef*L2D)) - (cWg)/(t_l))-c_f;

