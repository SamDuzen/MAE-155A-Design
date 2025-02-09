clear all; close all; clc;

%% Previous Parameter Input
Wg = 65000; %gross weight [N]
WL = 3640; %wing loading [N/m^2]
Ma = 1.6; %design mach number

%% New Set Aerodynamic Parameters
TaR = 0.24; %taper ratio (based on F22 and F35)
ThR = 0.045; %thickness ratio (based on historical trends from textbook)

a = 5.416; %historical aspect ratio parameter
c = -0.622; %historical aspect ratio parameter

WTA = -3; %wing twist angle [deg]
WIA = 0; %wing incidence angle [deg]
DA = 0; %dihedral angle [deg]


%% Calculations
AR = a*Ma^c; %aspect ratio
Sw_LE = 90 - asind(1/Ma); %leading edge sweep angle
Sw_QC = atand(tand(Sw_LE) - (1-TaR)/(AR*(1+TaR))); %quarter chord sweep angle
S = Wg/WL; %reference wing area [m^2]
b = sqrt(AR*S); %span [m]
C_Root = (2*S)/(b*(1+TaR)); %root chord [m]
C_Tip = TaR*C_Root; %tip chord [m]
MAC = (2/3)*((1+TaR+TaR^2)/(1+TaR))*C_Root;
Y_MAC = (b/6)*((1+2*TaR)/(1+TaR));

%% Table
Parameter_Name = {'Reference wing area [m^2]';'Span [m]';'Root chord [m]'; 'Tip chord [m]'; 'Mean aerodynamic chord'; 'MAC spanwise location';'Leading edge sweep angle [deg]'; 'Quarter chord sweep angle [deg]'; 'Taper ratio'; 'Thickness ratio'; 'Aspect ratio'; 'Twist angle [deg]'; 'Wing incidence angle [deg]'; 'Dihedral angle [deg]'};
Parameter_Symbol = {'S';'b';'C_root'; 'C_tip'; 'MAC'; 'Y'; 'Lambda_LE'; 'Lambda_c/4'; 'lambda'; 't/c'; 'AR'; 'theta'; 'i'; 'gamma'};
Parameter_Value = [S;b;C_Root;C_Tip;MAC;Y_MAC;Sw_LE;Sw_QC;TaR;ThR;AR;WTA;WIA;DA];

Column_Name = {'Parameter Name'; 'Symbol'; 'Value'};

Wing_Table = table(Parameter_Name, Parameter_Symbol, Parameter_Value,'VariableNames',Column_Name)

%% Tail Sizing





