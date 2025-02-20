clear all; close all; clc;

%% Adjustable Parameters

%General Sizing
Wg = 122917.9; %gross weight [N]
WL = 3640; %wing loading [N/m^2]
Ma = 1.6; %design mach number
L_Fus = 16.5; %fuselage length [m]

%Main Wing Aerodynamic Parameters
TaR = 0.24; %taper ratio (based on F22 and F35)
ThR = 0.045; %thickness ratio (based on historical trends from textbook)
a = 5.416; %historical aspect ratio parameter
c = -0.622; %historical aspect ratio parameter
WTA = -3; %wing twist angle [deg]
WIA = 0; %wing incidence angle [deg]
DA = 0; %dihedral angle [deg]

%Tail Parameters
L_VT = 0.5*L_Fus; %Tail moment arm as a percent of fuselage length
L_HT = 0.5*L_Fus; %Tail moment arm as a percent of fuselage length


%% Main Wing Sizing Calculations
AR = a*Ma^c; %aspect ratio
Sw_LE = 90 - asind(1/Ma); %leading edge sweep angle
Sw_QC = atand(tand(Sw_LE) - (1-TaR)/(AR*(1+TaR))); %quarter chord sweep angle
S = Wg/WL; %reference wing area [m^2]
b = sqrt(AR*S); %span [m]
C_Root = (2*S)/(b*(1+TaR)); %root chord [m]
C_Tip = TaR*C_Root; %tip chord [m]
MAC = (2/3)*((1+TaR+TaR^2)/(1+TaR))*C_Root;
Y_MAC = (b/6)*((1+2*TaR)/(1+TaR));

%% Tail Sizing

%Tail Parameters (Adjustable)
C_HT = 0.4; %Horizontal tail volume coefficient (jet fighter) - Historical
C_VT = 0.12; %Vertical tail volume coefficient (jet fighter) - Historical
AR_HT = 2; %Horizontal Tail aspect ratio (historical)
TaR_HT = 0.24; %Horizontal Tail taper ratio (historical)
AR_VT = 2; %Vertical Tail aspect ratio (historical)
TaR_VT = 0.24; %Vertical Tail taper ratio (historical)

%Tail Area
S_VT = (C_VT*b*S)/(L_VT); %Vertical Tail Area [m^2}
S_HT = (C_HT*MAC*S)/(L_HT); %Horizontal Tail Area [m^2]

%Horizontal Tail Calculations
b_HT = sqrt(AR_HT*S_HT); %span [m]
C_Root_HT = (2*S_HT)/(b_HT*(1+TaR_HT)); %root chord [m]
C_Tip_HT = TaR_HT*C_Root_HT; %tip chord [m]
MAC_HT = (2/3)*((1+TaR_HT+TaR_HT^2)/(1+TaR_HT))*C_Root_HT;
Y_MAC_HT = (b_HT/6)*((1+2*TaR_HT)/(1+TaR_HT));

%Vertical Tail Calculations
b_VT = sqrt(AR_VT*S_VT); %span [m]
C_Root_VT = (2*S_VT)/(b_VT*(1+TaR_VT)); %root chord [m]
C_Tip_VT = TaR_VT*C_Root_VT; %tip chord [m]
MAC_VT = (2/3)*((1+TaR_VT+TaR_VT^2)/(1+TaR_VT))*C_Root_VT;
Y_MAC_VT = (b_VT/6)*((1+2*TaR_VT)/(1+TaR_VT));


%% Table
Parameter_Name = {'Reference wing area [m^2]';'Span [m]';'Root chord [m]'; 'Tip chord [m]'; 'Mean aerodynamic chord'; 'MAC spanwise location';'Leading edge sweep angle [deg]'; 'Quarter chord sweep angle [deg]'; 'Taper ratio'; 'Thickness ratio'; 'Aspect ratio'; 'Twist angle [deg]'; 'Wing incidence angle [deg]'; 'Dihedral angle [deg]';'Vertical Tail Area [m^2]';'Horizontal Tail Area [m^2]'; 'Horizontal Tail Span [m]'; 'Horizontal Tail Root Chord [m]'; 'Horizontal Tail Tip Chord [m]'; 'Horizontal Tail Aspect Ratio'; 'Horizontal Tail Taper Ratio'; 'MAC (Horizontal Tail)'; 'MAC (Horizontal Tail) Spanwise Location';'MAC (Vertical Tail)'; 'MAC (Vertical Tail) Spanwise Location';'Vertical Tail Span [m]'; 'Vertical Tail Root Chord [m]'; 'Vertical Tail Tip Chord [m]'; 'Vertical Tail Aspect Ratio'; 'Vertical Tail Taper Ratio'};
Parameter_Symbol = {'S';'b';'C_root'; 'C_tip'; 'MAC'; 'Y'; 'Lambda_LE'; 'Lambda_c/4'; 'lambda'; 't/c'; 'AR'; 'theta'; 'i'; 'gamma';'S_VT';'S_HT';'b_H_tail'; 'C_root,Htail'; 'C_tip,Htail';'AR_HT';'lambda_HT';'MAC_HT';'Y_HT';'MAC_VT';'Y_VT';'b_V_tail'; 'C_root,Vtail'; 'C_tip,Vtail';'AR_VT';'lambda_VT'};
Parameter_Value = [S;b;C_Root;C_Tip;MAC;Y_MAC;Sw_LE;Sw_QC;TaR;ThR;AR;WTA;WIA;DA;S_VT;S_HT;b_HT;C_Root_HT;C_Tip_HT;AR_HT;TaR_HT;MAC_HT;Y_MAC_HT;MAC_VT;Y_MAC_VT;b_VT;C_Root_VT;C_Tip_VT;AR_VT;TaR_VT];

Column_Name = {'Parameter Name'; 'Symbol'; 'Value'};

Wing_Table = table(Parameter_Name, Parameter_Symbol, Parameter_Value,'VariableNames',Column_Name)

%% Export Data
writematrix(Parameter_Value,'Wing_Parameters.csv')