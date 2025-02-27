function [CG,NP_Sub,NP_Super,SM_Sub,SM_Super] = Weight_Balance(Wg,Ma_cruise)

%% Global Parameters
Ma = 1.6; %Design Mach Number

%% Pull Specified Parameters
t = readtable("Parameters.tsv", "FileType","text",'Delimiter', '\t');
    dCL_MW = t.(2)(17);
    dCL_HT = t.(2)(18);
    dCL_VT = t.(2)(19);
    X_MW = t.(2)(20); %Leading edge of main wing
    X_HT = t.(2)(21); %Leading edge of horizontal tail wing
    X_VT = t.(2)(22); %Leading edge of vertical tail wing
    CG_MW = t.(2)(23); %Main wing CG
    CG_HW = t.(2)(24); %Horizontal wing CG
    CG_VW = t.(2)(25); %Vertical wing CG
    CG_FU = t.(2)(26); %Fuselage CG
    CG_EN = t.(2)(27); %Engine CG
    CG_LG = t.(2)(28); %Landing Gear CG
    FL = t.(2)(29); %specified fuselage length [m]
    Sf_Wet = t.(2)(30); %fuselage wetted area [m^2] - CURRENTLY THE F-35 AS AN ESTIMATE
    W_eng_o = t.(2)(31)*t.(2)(2); %dry weight of engine [N]

%% Read Updated Wing Sizing Parameters
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
   

%% Approximate Empty Weight Buildup

%Approximate table values
W_Wing = 44*S; %Wing mass [kg]
W_Tail_H = 20*S_HT; %Horizontal tail mass [kg]
W_Tail_V = 20*S_VT; %Vertical tail mass [kg]
W_Fus = 23*Sf_Wet; %Fuselage weight [kg] - FIND Sf_Wet

W_LG = 0.033*Wg; %Landing gear weight [N]
W_AEE = 0.17*Wg; %"All-else empty" weight [N]

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


end