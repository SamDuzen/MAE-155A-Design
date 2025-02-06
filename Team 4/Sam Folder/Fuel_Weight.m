function [outputArg1,outputArg2] = Fuel_Weight(Range1,Range2,time1,time2,LD)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
Warmup_Fraction = 0.970;
Climb_Fraction = 0.985;
Landing_Fraction = 0.995;

C_cruise = 25.5; %mg/Ns
C_loiter = 22.7; %mg/Ns

%% Calculations

R1 = (1852)*Range1; %meters (from nautical miles)
R2 = (1852)*Range2; %meters (from nautical miles)
t1 = (3600)*time1; %seconds (from hours)
t2 = (3600)*time2; %seconds (from hours)


WR_1 = Warmup_Fraction; %Warmup and takeoff
WR_2 = Climb_Fraction;
WR_3 = exp((-R1*C_cruise)/(V*LD));
WR_4 = exp(-(t1*C_loiter)/(LD));
WR_5 = 0;
WR_6 = 0;
WR_7 = Landing_Fraction; %Landing
WR_8 = 0;

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end