clear all; close all; clc;

%% Adjustable Parameters
%General sizing
Wg = 65000; %gross weight [N]
T2W = 0.56; %thrust to weight ratio
WL = 3640; %wing loading [N/m^2]

%General mission
Ma = 1.6; %design mach number
Ma_cruise = 0.85;
e = 0.7; %efficiency factor
AR = 4.0431; %aspect ratio
n_missile = 4; %number of missiles

%General propulsion
n_eng = 1; %number of engines being used

%General Aerodynamics
Cdo_c = 0.0468; %parasitic drag at cruise

%% Improved Payload Weight Estimate

%Government Furnished Equipment (excluding weapons)
W_gov = 100+10+100+50+50+450+300+100;
if n_eng == 1
    W_gov = W_gov - 80;
end

%Weapons
W_missiles = 327*n_missile; %total missile weight [lb]
W_cannon = 275+300+500*(0.58); %loaded cannon weight [lb]

%Additional Payloads
    %Nothing for now

Wp = W_gov+W_missiles+W_cannon;
    Wp_BU = Wp; %storing british units weight
    Wp = Wp*4.4482; %conversion: [lb] --> [N]

%% Improved Empty Weight Fraction

%Coefficients (for a jet fighter)
a = 0;
b = 4.28;
C1 = -0.1; 
C2 = 0.1;
C3 = 0.2;
C4 = -0.24;
C5 = 0.11;

%New empty weight fraction estimate
EWF = a + b*(Wg^C1)*(AR^C2)*(T2W^C3)*(WL^C4)*(Ma^C5);

%% Improved L/D Estimate (For cruise/loiter)

%Dynamic Pressure @ Cruise - 9500 m (31000 ft)
den_95 = 0.43966; %[kg/m^3]
a_95 = 301.7; %[m/s]
v_95 = Ma_cruise*a_95; %cruise velocity
q = 0.5*den_95*v_95^2; %dynamic pressure [N/m^2]

%L/D Estimation
L2D = 1/(((q*Cdo_c)/(WL))+WL*(1/(q*pi*AR*e)));


%% Improved Fuel Weight Estimate

%% Fuselage

%Estimation Parameters (Jet Fighter)
    a = 0.389;
    c = 0.39;
FL = a*Wg^c;