clear all; close all; clc;

%% Parameters
payload = 5000; %lb
max_n = 7; %Maximum design load factor
AR = 2.5;

%Empty Weight Fraction
a = 1.47; %UCAV
c = -0.16; %UCAV
Kvs = 1; %Fixed Sweep (change to 1.04 for variable sweep)

%% Corner Speed Constraint

% Given Parameters
TR_max = 18; %deg/s
TRN_M1 = 1.2; %Mach Number

% Natural Parameters
a_35 = 659; %speed of sound (mph) at 35000ft
g = 9.81; %m/s^2
density_35 = 0.000736; %air density (standard atmosphere) @ 35000 ft

% Calculation
a_35 = a_35/2.237; %conversion from mph to m/s
    V_35 = TRN_M1*a_35; %Mach 1.2 velocity @ 35000 ft
TR_max_rad = (TR_max/360)*(2*pi); %conversion from deg/s to rad/s
n_turn = sqrt(((TR_max_rad*V_35)/(g))^2 + 1);

% Maximum combat wing loading

WL_Combat = ((0.5*density_35*V_35^2)*0.6)/(n_turn);









%% Functions