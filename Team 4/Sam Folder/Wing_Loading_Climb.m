clear all; close all; clc;

%% Parameters
Ma = 0.9;

%% Performance Requirements

%1-g Specific Excess Power - Military Thrust
n1_MiT_SL = 200; %Sea Level - [ft/sec]
n1_MiT_15 = 50; %15000 ft - [ft/sec]

%1-g Specific Excess Power - Maximum Thrust
n1_MaT_SL = 700; %Sea Level - [ft/sec]
n1_MaT_15 = 400; %15000 ft - [ft/sec]

%5-g Specific Excess Power - Maximum Thrust
n5_MaT_SL = 300; %Sea Level - [ft/sec]
n5_MaT_15 = 50; %15000 ft - [ft/sec]

%% Altitude Conditions

%Sea Level
den_SL = 0.00238; %density [slugs/ft^3]
a_SL = 1116.5; %speed of sound [ft/s]

%15000 ft
den_15 = 0.00150; %density [slugs/ft^3]
a_15 = 1057.4; %speed of sound [ft/s]

%% Velocity Calculations

%Sea Level
V_SL = a_SL*Ma;

%15000 ft
V_15 = a_15*Ma;

%% Climb Angle Calculations

%1-g Climb Angle - Military Thrust
CA(1) = asind(n1_MiT_SL/V_SL); %Sea Level - [deg]
CA(2) = asind(n1_MiT_15/V_15); %15000 ft - [deg]

%1-g Climb Angle - Maximum Thrust
CA(3) = asind(n1_MaT_SL/V_SL); %Sea Level - [deg]
CA(4) = asind(n1_MaT_15/V_15); %15000 ft - [deg]

%5-g Climb Angle - Maximum Thrust
CA(5) = asind(n5_MaT_SL/V_SL); %Sea Level - [deg]
CA(6) = asind(n5_MaT_15/V_15); %15000 ft - [deg]






