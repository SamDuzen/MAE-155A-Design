clear all; close all; clc;

%% Adjustable Parameters
%General sizing
Gross_Weight = 60000; %gross weight [N]
T2W = 0.56; %thrust to weight ratio
    T2W_dash = T2W; %change if needed
WL = 3640; %wing loading [N/m^2]
    WL_Cruise = WL; %change if needed
    WL_dash = WL; %change if needed

%General mission
Ma = 1.6; %design mach number
Ma_cruise = 0.85;
e = 0.7; %efficiency factor
    e_combat = e; %change if needed
g = 9.81;
AR = 4.0431; %aspect ratio
n_missile = 4; %number of missiles

%General propulsion
n_eng = 2; %number of engines being used
TSFC = 0.71; %Thrust specific fuel consumption (british units)
    TSFC_cruise = TSFC; %change if needed
    TSFC_dash = TSFC; %change if needed
W_eng_o = 9804*2; %dry weight of engine [N]

%General Aerodynamics
Cdo_c = 0.0468; %parasitic drag at cruise
Cdo_d = 0.0493; %parasitic drag at dash speed
dCL_MW = (0.7-0.2)/5*360/(2*pi); %Lift coefficient slope for main wing
dCL_HT = 0.1*360/(2*pi); %Lift coefficient slope for horizontal tail
dCL_VT = 0.075*360/(2*pi); %Lift coefficient slope for vertical tail

%Fuselage
FL_s = 16.5; %specified fuselage length [m]
    FusF = 1; %fuselage factor (set to 1 if fuselage length is specified, set to 0 if it should be calculated)
Sf_Wet = 36.1; %fuselage wetted area [m^2] - CURRENTLY THE F-35 AS AN ESTIMATE

%Positioning (from front of Fuselage) - [%]
X_MW = 0.25; %Leading edge of main wing
X_HT = 0.8; %Leading edge of horizontal tail wing
X_VT = 0.8; %Leading edge of vertical tail wing
CG_MW = 0.4471; %Main wing CG
CG_HW = 0.8832; %Horizontal wing CG
CG_VW = 0.8864; %Vertical wing CG
CG_FU = 0.55; %Fuselage CG
CG_EN = 0.7; %Engine CG
CG_LG = 0.65; %Landing Gear CG

%% Import Wing Data
Wing_Parameters = csvread("Wing_Parameters.csv");
    S = Wing_Parameters(1);
    b = Wing_Parameters(2);
    C_Root = Wing_Parameters(3);
    C_Tip = Wing_Parameters(4);
    MAC = Wing_Parameters(5);
    Y_MAC = Wing_Parameters(6);
    Sw_LE = Wing_Parameters(7);
    Sw_QC = Wing_Parameters(8);
    TaR = Wing_Parameters(9);
    ThR = Wing_Parameters(10);
    AR = Wing_Parameters(11);
    WTA = Wing_Parameters(12);
    WIA = Wing_Parameters(13);
    DA = Wing_Parameters(14);
    S_VT = Wing_Parameters(15);
    S_HT = Wing_Parameters(16);
    b_HT = Wing_Parameters(17);
    C_Root_HT = Wing_Parameters(18);
    C_Tip_HT = Wing_Parameters(19);
    AR_HT = Wing_Parameters(20);
    TaR_HT = Wing_Parameters(21);
    MAC_HT = Wing_Parameters(22);
    Y_MAC_HT = Wing_Parameters(23);
    MAC_VT = Wing_Parameters(24);
    Y_MAC_VT = Wing_Parameters(25);

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

%% Improved L/D Estimate (For cruise/loiter)

%Dynamic Pressure @ Cruise - 9500 m (31000 ft)
den_95 = 0.43966; %[kg/m^3]
a_95 = 301.7; %[m/s]
v_cruise = Ma_cruise*a_95; %cruise velocity
q_cruise = 0.5*den_95*v_cruise^2; %dynamic pressure [N/m^2]

%L/D Estimation
L2D_cruise = 1/(((q_cruise*Cdo_c)/(WL_Cruise))+WL_Cruise*(1/(q_cruise*pi*AR*e)));

%% L/D Estimate (For dash)
    %no known better method yet so for now...
a_dash = 297.4; %speed of sound at 10.5km (about 35000 ft) [m/s]
v_dash = Ma*a_dash;
den_dash = 0.388857; %[kg/m^3]
q_dash = 0.5*den_dash*v_dash^2;

L2D_dash = 1/(((q_dash*Cdo_d)/(WL_dash))+WL_dash*(1/(q_dash*pi*AR*e)));

%% L/D Estimate (for combat)
n_max = 7; %Instantaneous turn rate
L2D_combat = 1/(((q_dash*Cdo_d)/(n_max*WL_dash))+WL_dash*n_max*(1/(q_dash*pi*AR*e_combat)));

%% Start Numerical Solver
Wg(1) = Gross_Weight;
error(1) = 20000;
i = 1;
while abs(error(i)) > 10

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
EWF = a + b*(Wg(i)^C1)*(AR^C2)*(T2W^C3)*(WL^C4)*(Ma^C5);

%% Improved Fuel Weight Estimate
W(1) = Wg(i); %initial weight

%Engine Start/Taxi/Takeoff
Wfr(1) = 0.98; %historical estimate
W(2) = Wfr(1)*W(1); %historical estimate
Wfu(1) = (1-Wfr(1))*W(1); %fuel weight

%Climb and Accelerate
Wfr(2) = 1.0065 - 0.0325*Ma_cruise; %accelerate up to cruise
W(3) = Wfr(2)*W(2);
Wfu(2) = (1-Wfr(2))*W(2);

%Cruise - 300nm
R1 = 300*1852; %range in meters
if i==1
    TSFC_cruise = TSFC_cruise/3600/2.205/4.448;
end
Wfr(3) = exp(-(TSFC_cruise*g*R1)/(v_cruise*L2D_cruise));
W(4) = Wfr(3)*W(3);
Wfu(3) = (1-Wfr(3))*W(3);

%Loiter - 4hrs
t1 = 4*60*60; %loiter time [sec]
Wfr(4) = exp(-(TSFC_cruise*g*t1)/(L2D_cruise));
W(5) = Wfr(4)*W(4);
Wfu(4) = (1-Wfr(4))*W(4);

%Accelerate up to max speed
Wfr(5) = (0.991-0.007*Ma-0.01*Ma^2)/(1.0065-0.0325*Ma_cruise);
W(6) = Wfr(5)*W(5);
Wfu(5) = (1-Wfr(5))*W(5);

%Dash - 100nm
if i==1
    TSFC_dash = TSFC_dash/3600/2.205/4.448;
end
R2 = 100*1852; %range in meters
Wfr(6) = exp(-(TSFC_dash*g*R2)/(v_dash*L2D_dash));
W(7) = Wfr(6)*W(6);
Wfu(6) = (1-Wfr(6))*W(6);

%Combat
v_turn_1 = 1.2*a_dash; %turn speed for turn 1
v_turn_2 = 0.9*a_dash; %turn speed for turn 2
d1 = (2*pi*v_turn_1)/(g*sqrt(n_max^2-1)); %combat duration for turn 1
d2 = (2*pi*v_turn_2)/(g*sqrt(n_max^2-1)); %combat duration for turn 2
d = d1+d2; %total combat duration
Wfr(7) = 1 - TSFC_dash*g*T2W_dash*d;
W(8) = Wfr(7)*W(7);
    W(8) = W(8) - W_missiles; %all missiles fired
Wfu(7) = (1-Wfr(7))*W(7);

%Cruise - 400nm
R2 = 400*1852; %range in meters
Wfr(8) = exp(-(TSFC_cruise*g*R2)/(v_cruise*L2D_cruise));
W(9) = Wfr(8)*W(8);
Wfu(8) = (1-Wfr(8))*W(8);

%Descent
Wfr(9) = 0.9925; %historical estimation
W(10) = Wfr(9)*W(9);
Wfu(9) = (1-Wfr(9))*W(9);

%Landing/Taxi
Wfr(10) = 0.995; %historical estimation
W(11) = Wfr(10)*W(10);
Wfu(10) = (1-Wfr(10))*W(10);

%Reserves: 30 minutes loiter
t2 = 0.5*60*60; %loiter time [sec]
Wfr(11) = exp(-(TSFC_cruise*g*t2)/(L2D_cruise));
W(12) = Wfr(11)*W(11);
Wfu(11) = (1-Wfr(11))*W(11);

%Total Fuel Weight
Wf = sum(Wfu);
Wf = 1.06*Wf; %additional allowance for unusable fuel


%% Gross Weight
Wg_calc(i) = Wp+Wf+EWF*Wg(i);

%% End Numerical Solver
error(i+1) = abs(Wg_calc(i)-Wg(i));
if error(i+1) < error(i)
    Wg(i+1) = Wg(i)+0.05*abs(error(i+1));
elseif error(i+1) > error(i)
    Wg(i+1) = Wg(i)-0.05*abs(error(i+1));
end

i = i+1;
end

Wg_calc(end) = [];
fprintf('The gross weight converges at %d Newtons\n',Wg_calc(end))

%% Fuel Volume
FD = 6.7; %fuel density [lb/gal]
FD = FD*119.8266; %conversion: [lb/gal] --> [kg/m^3]
FV = (Wf/g)/FD; %fuel volume [m^3]
fprintf('The fuel volume is %d m^3\n',FV)

%% Payload Volume

%Government Furnished Equipment (excluding weapons)
V_gov = 3+0.5+3+1+2+6+4+2; %[ft^3]
if n_eng == 1
    V_gov = V_gov - 1;
end

%Weapons
V_missiles = pi*(0.6/2)^2*12*n_missile; %total missile volume [ft^3]
V_cannon_int = pi*(25/12/2)^2*(25/12); %internal cannon volume [ft^3]
V_cannon_ext = pi*(10/12/2)^2*(74/12); %external cannon volume [ft^3]

%Additional Payloads
    %Nothing for now

Vp_int = V_gov+V_cannon_int; %total payload internal volume [ft^3]
    Vp_int = Vp_int/35.315; %conversion: [ft^3] --> [m^3]
Vp_ext = V_missiles+V_cannon_ext; %total payload external volume [ft^3]
    Vp_ext = Vp_ext/35.315; %conversion: [ft^3] --> [m^3]

%Total Volume
fprintf('The total used internal volume is %d m^3\n',FV+Vp_int)
fprintf('The total used external volume is %d m^3\n',Vp_ext)

%% Fuselage

%Estimation Parameters (Jet Fighter)
    a = 0.389;
    c = 0.39;
FL = a*Wg_calc(end)^c; 

if FusF == 1
    FL = FL_s; %set to specified length
end

%% Approximate Empty Weight Buildup

%Approximate table values
W_Wing = 44*S; %Wing mass [kg]
W_Tail_H = 20*S_HT; %Horizontal tail mass [kg]
W_Tail_V = 20*S_VT; %Vertical tail mass [kg]
W_Fus = 23*Sf_Wet; %Fuselage weight [kg] - FIND Sf_Wet

W_LG = 0.033*Wg_calc(end); %Landing gear weight [N]
W_AEE = 0.17*Wg_calc(end); %"All-else empty" weight [N]

W_eng_i = 1.3*W_eng_o; %Installed engine weight [N]

%Total weight check
W_tot = (W_Wing+W_Tail_H+W_Tail_V+W_Fus)*9.81 + W_LG + W_AEE + W_eng_i;

%% Center of Gravity

X_CG = (9.81*(W_Wing*CG_MW + W_Tail_H*CG_HW + W_Tail_V*CG_VW + W_Fus*CG_FU)+W_LG*CG_LG+W_eng_i*CG_EN+W_AEE*0.5)/(W_tot);

%% Aerodynamic Center

%Airspeed correction factor
dX_AC_Sub = 0.26*(Ma_cruise-0.4)^2.5; %0.4<Ma<1.1
dX_AC_Super = 0.112 - 0.004*Ma; %M>1.1

%Main Wing
X_MAC_QC_MW = X_MW*FL + Y_MAC*tand(Sw_LE) + MAC/4; %Quarter chord of MAC from front of aircraft
    X_AC_MW_Sub = X_MAC_QC_MW + dX_AC_Sub*sqrt(S); %Subsonic
    X_AC_MW_Super = X_MAC_QC_MW + dX_AC_Super*sqrt(S); %Supersonic

%Horizontal Tail Wing
X_MAC_QC_HT = X_HT*FL + Y_MAC_HT*tand(Sw_LE) + MAC_HT/4; %Quarter chord of MAC from front of aircraft
    X_AC_HT_Sub = X_MAC_QC_HT + dX_AC_Sub*sqrt(S_HT); %Subsonic
    X_AC_HT_Super = X_MAC_QC_HT + dX_AC_Super*sqrt(S_HT); %Supersonic

%Verticle Tail Wing
X_MAC_QC_VT = X_VT*FL + Y_MAC_VT*tand(Sw_LE) + MAC_VT/4; %Quarter chord of MAC from front of aircraft
    X_AC_VT_Sub = X_MAC_QC_VT + dX_AC_Sub*sqrt(S_VT); %Subsonic
    X_AC_VT_Super = X_MAC_QC_VT + dX_AC_Super*sqrt(S_VT); %Supersonic

%% Neutral Point

%Subsonic
X_NP_Sub = (dCL_MW*S*X_AC_MW_Sub + dCL_HT*S_HT*X_AC_HT_Sub + dCL_VT*S_VT*X_AC_VT_Sub)/(dCL_MW*S + dCL_HT*S_HT + dCL_VT*S_VT);

%Supersonic
X_NP_Super = (dCL_MW*S*X_AC_MW_Super + dCL_HT*S_HT*X_AC_HT_Super + dCL_VT*S_VT*X_AC_VT_Super)/(dCL_MW*S + dCL_HT*S_HT + dCL_VT*S_VT);

%% Converting to % of MAC from Leading Edge

%MAC Leading Edge Locations
X_MAC_MW = X_MW*FL + Y_MAC*tand(Sw_LE);
X_MAC_HT = X_HT*FL + Y_MAC_HT*tand(Sw_LE);
X_MAC_VT = X_VT*FL + Y_MAC_VT*tand(Sw_LE);

%Center of Gravity
CG = abs(X_CG*FL-X_MAC_MW)/MAC;

%Neutral Point - Subsonic
NP_Sub = abs(X_NP_Sub-X_MAC_MW)/MAC;

%Neutral Point - Supersonic
NP_Super = abs(X_NP_Super-X_MAC_MW)/MAC;

%% Static Margin

%Subsonic
SM_Sub = (X_NP_Sub-X_CG*FL)/MAC

%Supersonic
SM_Super = (X_NP_Super-X_CG*FL)/MAC




%% Start Numerical Solver [NOTE: THIS IS JUST TO TEST. EVERYTHING BELOW CAN BE DELTED]
Wgt(1) = Gross_Weight;
error(1) = 20000;
i = 1;
while abs(error(i)) > 10

%% Wing Sizing

%Main Wing
St(i) = Wgt(i)/WL;
bt(i) = sqrt(AR*St(i)); %span [m]
C_Roott(i) = (2*St(i))/(bt(i)*(1+TaR)); %root chord [m]
MACt(i) = (2/3)*((1+TaR+TaR^2)/(1+TaR))*C_Roott(i);

%Horizontal Tail
C_HT = 0.4; %Horizontal tail volume coefficient (jet fighter) - Historical
L_HT = 0.5*FL; %Tail moment arm as a percent of fuselage length
S_HTt(i) = (C_HT*MACt(i)*St(i))/(L_HT); %Horizontal Tail Area [m^2]

%Vertical Tail
C_VT = 0.12; %Vertical tail volume coefficient (jet fighter) - Historical
L_VT = 0.5*FL; %Tail moment arm as a percent of fuselage length
S_VTt(i) = (C_VT*bt(i)*St(i))/(L_VT); %Vertical Tail Area [m^2}
    
    %% Approximate Empty Weight Buildup

%Approximate table values
W_Wing = 44*St(i); %Wing mass [kg]
W_Tail_H = 20*S_HTt(i); %Horizontal tail mass [kg]
W_Tail_V = 20*S_VTt(i); %Vertical tail mass [kg]
W_Fus = 23*Sf_Wet; %Fuselage weight [kg] - FIND Sf_Wet

W_LG = 0.033*Wgt(i); %Landing gear weight [N]
W_AEE = 0.17*Wgt(i); %"All-else empty" weight [N]

W_eng_i = 1.3*W_eng_o; %Installed engine weight [N]

%Total weight check
W_tot_t(i) = (W_Wing+W_Tail_H+W_Tail_V+W_Fus)*9.81 + W_LG + W_AEE + W_eng_i;

%% Improved Empty Weight Fraction

%Coefficients (for a jet fighter)
EW = W_tot_t(i);

%% Improved Fuel Weight Estimate
W(1) = Wgt(i); %initial weight

%Engine Start/Taxi/Takeoff
Wfr(1) = 0.98; %historical estimate
W(2) = Wfr(1)*W(1); %historical estimate
Wfu(1) = (1-Wfr(1))*W(1); %fuel weight

%Climb and Accelerate
Wfr(2) = 1.0065 - 0.0325*Ma_cruise; %accelerate up to cruise
W(3) = Wfr(2)*W(2);
Wfu(2) = (1-Wfr(2))*W(2);

%Cruise - 300nm
R1 = 300*1852; %range in meters
Wfr(3) = exp(-(TSFC_cruise*g*R1)/(v_cruise*L2D_cruise));
W(4) = Wfr(3)*W(3);
Wfu(3) = (1-Wfr(3))*W(3);

%Loiter - 4hrs
t1 = 4*60*60; %loiter time [sec]
Wfr(4) = exp(-(TSFC_cruise*g*t1)/(L2D_cruise));
W(5) = Wfr(4)*W(4);
Wfu(4) = (1-Wfr(4))*W(4);

%Accelerate up to max speed
Wfr(5) = (0.991-0.007*Ma-0.01*Ma^2)/(1.0065-0.0325*Ma_cruise);
W(6) = Wfr(5)*W(5);
Wfu(5) = (1-Wfr(5))*W(5);

%Dash - 100nm
R2 = 100*1852; %range in meters
Wfr(6) = exp(-(TSFC_dash*g*R2)/(v_dash*L2D_dash));
W(7) = Wfr(6)*W(6);
Wfu(6) = (1-Wfr(6))*W(6);

%Combat
v_turn_1 = 1.2*a_dash; %turn speed for turn 1
v_turn_2 = 0.9*a_dash; %turn speed for turn 2
d1 = (2*pi*v_turn_1)/(g*sqrt(n_max^2-1)); %combat duration for turn 1
d2 = (2*pi*v_turn_2)/(g*sqrt(n_max^2-1)); %combat duration for turn 2
d = d1+d2; %total combat duration
Wfr(7) = 1 - TSFC_dash*g*T2W_dash*d;
W(8) = Wfr(7)*W(7);
    W(8) = W(8) - W_missiles; %all missiles fired
Wfu(7) = (1-Wfr(7))*W(7);

%Cruise - 400nm
R2 = 400*1852; %range in meters
Wfr(8) = exp(-(TSFC_cruise*g*R2)/(v_cruise*L2D_cruise));
W(9) = Wfr(8)*W(8);
Wfu(8) = (1-Wfr(8))*W(8);

%Descent
Wfr(9) = 0.9925; %historical estimation
W(10) = Wfr(9)*W(9);
Wfu(9) = (1-Wfr(9))*W(9);

%Landing/Taxi
Wfr(10) = 0.995; %historical estimation
W(11) = Wfr(10)*W(10);
Wfu(10) = (1-Wfr(10))*W(10);

%Reserves: 30 minutes loiter
t2 = 0.5*60*60; %loiter time [sec]
Wfr(11) = exp(-(TSFC_cruise*g*t2)/(L2D_cruise));
W(12) = Wfr(11)*W(11);
Wfu(11) = (1-Wfr(11))*W(11);

%Total Fuel Weight
Wf = sum(Wfu);
Wf = 1.06*Wf; %additional allowance for unusable fuel


%% Gross Weight
Wg_calct(i) = Wp+Wf+EW;

%% End Numerical Solver
error(i+1) = abs(Wg_calct(i)-Wgt(i));
if error(i+1) < error(i)
    Wgt(i+1) = Wgt(i)-0.05*abs(error(i+1));
elseif error(i+1) > error(i)
    Wgt(i+1) = Wgt(i)+0.05*abs(error(i+1));
end

i = i+1;
end

Wg_calct(end) = [];
fprintf('The gross weight converges at %d Newtons\n',Wg_calct(end))

