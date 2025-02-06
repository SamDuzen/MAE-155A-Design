clear all; close all; clc;

%% Sizing Parameters

%time/range (mission segments)
R1 = 300; %nm
t2 = 4; %hours
R3 = 100; %nm
t4 = 0.5; %hours

R_test = 1600000; %meters

%engine
Ct = 0.00017;

%aircraft
LD = 10;
v = 1200; %mph
    v = v*0.44704;
a = 2.11; %coefficient (empty weight)
b = -0.13; %coefficient (empty weight)

%global
g = 9.81; %m/s

%payload
Wp = 2352; %lb
%% Fuel Weight

%Segment 1
R1 = R1*1852; %conversion to meters from nautical miles
WR_f1 = 1 - exp(-((R1*Ct*g)/(LD*v)));

%Segment 2
t2 = t2*60*60;
WR_f2 = 1 - WR_f1 - (1 - WR_f1)*exp(-((t2*Ct*g)/(LD)));

%Segment 3
R3 = R3*1852;
WR_f3 = 1 - WR_f1 - WR_f2 - (1 - WR_f1 - WR_f2)*exp(-((R3*Ct*g)/(LD*v)));

%Segment 4
t4 = t4*60*60;
WR_f4 = 1 - WR_f1 - WR_f2 - WR_f3 - (1 - WR_f1 - WR_f2 - WR_f3)*exp(-((t4*Ct*g)/(LD)));

%Test
WR_ft = 1 - exp(-((R_test*Ct*g)/(LD*v)));

%% Payload Weight
Wp = Wp*0.4536; %kg

%% Empty Weight Fraction

