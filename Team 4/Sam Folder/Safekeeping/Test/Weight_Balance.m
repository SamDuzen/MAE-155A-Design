function [CG_WM,CG_NM,X_CG_WM,X_CG_NM,NP_Sub,NP_Super,SM_Sub_WM,SM_Sub_NM,SM_Super_WM,SM_Super_NM,d_display] = Weight_Balance(Wg,Ma_cruise,Wf,Wfu)

%% Global Parameters
Ma = 1.6; %Design Mach Number

%% Pull Specified Parameters
t = readtable("Parameters.tsv", "FileType","text",'Delimiter', '\t');
    n_missile = t.(2)(1);
    n_eng = t.(2)(2);
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
   

%% Payload Weight

%Payload locations (from the front of fuselage)
X_ICNIA = 185*0.0254;
X_DataBus = 345*0.0254;
X_INEWS = 165*0.0254;
X_VMS = 370*0.0254;
X_IRSTS = 65*0.0254;
X_AAR = 39.37*0.0254;
X_ES = 365*0.0254;
X_APU = 420*0.0254;

%Payload weights
W_ICNIA = 100*4.44822;
W_DataBus = 10*4.44822;
W_INEWS = 100*4.44822;
W_VMS = 50*4.44822;
W_IRSTS = 50*4.44822;
W_AAR = 450*4.44822;
if n_eng == 2
    W_ES = 300*4.44822;
else
    W_ES = 220*4.44822;
end
W_APU = 100*4.44822;

CG_Payload = (W_ICNIA*X_ICNIA + W_DataBus*X_DataBus + W_INEWS*X_INEWS + W_VMS*X_VMS + W_IRSTS*X_IRSTS + W_AAR*X_AAR + W_ES*X_ES + W_APU*X_APU)/(W_ICNIA + W_DataBus + W_INEWS + W_VMS + W_IRSTS + W_AAR + W_ES + W_APU);
CG_Payload = CG_Payload/FL; %Nondimensionalize
W_Payload = W_ICNIA + W_DataBus + W_INEWS + W_VMS + W_IRSTS + W_AAR + W_ES + W_APU;

%% Weapons
X_Gun = 0.4;
W_Gun = (300+275)*4.44822;

X_Missiles = CG_MW;
W_Missiles = 327*n_missile*4.44822;

CG_Weapons_WM = (W_Gun*X_Gun + W_Missiles*X_Missiles)/(W_Gun + W_Missiles); %With missiles
CG_Weapons_NM = X_Gun; %No missiles

W_Weapons_WM = W_Gun + W_Missiles; %With missiles
W_Weapons_NM = W_Gun; %No missiles
    
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

%Weight Table
Component_Name = {'Wing Weight';'Vertical Tail Weight'; 'Fuselage Weight'; 'Landing Gear Weight';'Engine';'Horizontal Tail Weight'; 'Payload'; 'Weapons'; 'All Else Empty'};
Component_Weight = [W_Wing*9.81; W_Tail_V*9.81; W_Fus*9.81; W_LG; W_eng_i; W_Tail_H*9.81; W_Payload; W_Weapons_WM; W_AEE];

Column_Name = {'Parameter Name'; 'Value'};

%Weight Charts
Weight_Table = table(Component_Name,Component_Weight,'VariableNames',Column_Name);

Weight_Table = table(Component_Name,Component_Weight);

p = piechart(Weight_Table,"Component_Weight","Component_Name");

% Modify font size for labels
for i = 1:length(p)
    if isprop(p(i), 'FontSize') % Check if it's a text object
        p(i).FontSize = 16;  % Adjust font size as needed
    end
end

p


%% Empty Weight Center of Gravity

CG_Empty = (9.81*(W_Wing*CG_MW + W_Tail_H*CG_HW + W_Tail_V*CG_VW + W_Fus*CG_FU)+W_LG*CG_LG+W_eng_i*CG_EN+W_AEE*0.5)/(W_tot);

X_CG_WM = (W_Payload*CG_Payload + W_Weapons_WM*CG_Weapons_WM + W_tot*CG_Empty)/(W_Payload + W_Weapons_WM + W_tot); %overall CG with missiles
X_CG_NM = (W_Payload*CG_Payload + W_Weapons_NM*CG_Weapons_NM + W_tot*CG_Empty)/(W_Payload + W_Weapons_NM + W_tot); %overall CG without missiles

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
CG_WM = abs(X_CG_WM*FL-X_MAC_MW)/MAC;
CG_NM = abs(X_CG_NM*FL-X_MAC_MW)/MAC;

%Neutral Point - Subsonic
NP_Sub = abs(X_NP_Sub-X_MAC_MW)/MAC;

%Neutral Point - Supersonic
NP_Super = abs(X_NP_Super-X_MAC_MW)/MAC;

%% Static Margin

%Subsonic
SM_Sub_WM = (X_NP_Sub-X_CG_WM*FL)/MAC; %with missiles
SM_Sub_NM = (X_NP_Sub-X_CG_NM*FL)/MAC; %without missiles

%Supersonic
SM_Super_WM = (X_NP_Super-X_CG_WM*FL)/MAC; %with missiles
SM_Super_NM = (X_NP_Super-X_CG_NM*FL)/MAC; %without missiles

%% Sensitivity Study (fuel weight)

Left_Bound_Input = 0.3; %[m]
    Left_Bound = Left_Bound_Input/FL; %Convert to percentage
Right_Bound_Input = 0.2; %[m]
    Right_Bound = Right_Bound_Input/FL; %Convert to percentage
Steps = 1000;

W_Overall_WM = W_Payload + W_Weapons_WM + W_tot;
W_Overall_NM = W_Payload + W_Weapons_NM + W_tot;
d(1) = -Left_Bound;
d_display(1) = -Left_Bound_Input;
placeholder1 = 0;
placeholder2 = 0;
for i = 1:Steps
    %Find Fuel CG
        CG_Fuel_WM = X_CG_WM + d(i);
        CG_Fuel_NM = X_CG_NM + d(i);
    %Find Overall CG
        CG_Overall_WM = (Wf*CG_Fuel_WM + W_Overall_WM*X_CG_WM)/(Wf + W_Overall_WM);
        CG_Overall_NM = (Wf*CG_Fuel_NM + W_Overall_NM*X_CG_NM)/(Wf + W_Overall_NM);
    %Find Static Margin
        SM_Fuel_Sub_WM(i) = (X_NP_Sub-CG_Overall_WM*FL)/MAC; %with missiles
        SM_Fuel_Sub_NM(i) = (X_NP_Sub-CG_Overall_NM*FL)/MAC; %without missiles
        SM_Fuel_Super_WM(i) = (X_NP_Super-CG_Overall_WM*FL)/MAC; %with missiles
        SM_Fuel_Super_NM(i) = (X_NP_Super-CG_Overall_NM*FL)/MAC; %with missiles
    %Report
        SM_Vector = [SM_Fuel_Sub_WM(i),SM_Fuel_Sub_NM(i),SM_Fuel_Super_WM(i),SM_Fuel_Super_NM(i)];
        if sum(SM_Vector >= 0.1) == 1
            Upper_Max = d_display(i);
        end
        if placeholder2 == 0
            if any(SM_Vector <= -0.1)
                Lower_Max = d_display(i);
                placeholder2 = placeholder2 + 1;
            end
        end
    %Update d
        d(i+1) = d(i) + (Left_Bound + Right_Bound)/Steps;
        d_display(i+1) = d_display(i) + (Left_Bound_Input + Right_Bound_Input)/Steps;
end

d(end) = [];
d_display(end) = [];

%Plot
figure()
hold on
plot(d_display,SM_Fuel_Sub_WM,'LineWidth',1)
plot(d_display,SM_Fuel_Sub_NM,'LineWidth',1)
plot(d_display,SM_Fuel_Super_WM,'LineWidth',1)
plot(d_display,SM_Fuel_Super_NM,'LineWidth',1)
yline(0.1,'--k')
yline(-0.1,'--k')
yline(0,'k')
xline(Lower_Max,'--r')
xline(Upper_Max,'--r')
ylabel('Static Margin')
xlabel('Fuel CG Distance [m]')
legend('Subsonic (w/ Missiles)','Subsonic (No Missiles)','Supersonic (w/ Missiles)','Supersonic (No Missiles)','Location','best')
title('100% Fuel Capacity')

%% Sensitivity Study (sloshing)
Steps = 1000;
Fuel_Step = 100;

W_Fuel(1) = (1/Fuel_Step)*Wf;
FW_Display(1) = (1/Fuel_Step)*100;
for j = 1:Fuel_Step

d_WM(1) = -X_CG_WM; %start at 0% FL
d_NM(1) = -X_CG_NM; %start at 0% FL
d_display_WM(1) = d_WM(1)*FL;
d_display_NM(1) = d_NM(1)*FL;
placeholder1 = 0;
placeholder2 = 0;
for i = 1:Steps
    %Find Fuel CG
        CG_Fuel_WM = X_CG_WM + d_WM(i);
        CG_Fuel_NM = X_CG_NM + d_NM(i);
    %Find Overall CG
        CG_Overall_WM = (W_Fuel(j)*CG_Fuel_WM + W_Overall_WM*X_CG_WM)/(W_Fuel(j) + W_Overall_WM);
        CG_Overall_NM = (W_Fuel(j)*CG_Fuel_NM + W_Overall_NM*X_CG_NM)/(W_Fuel(j) + W_Overall_NM);
    %Find Static Margin
        SM_Fuel_Sub_WM(i) = (X_NP_Sub-CG_Overall_WM*FL)/MAC; %with missiles
        SM_Fuel_Sub_NM(i) = (X_NP_Sub-CG_Overall_NM*FL)/MAC; %without missiles
        SM_Fuel_Super_WM(i) = (X_NP_Super-CG_Overall_WM*FL)/MAC; %with missiles
        SM_Fuel_Super_NM(i) = (X_NP_Super-CG_Overall_NM*FL)/MAC; %with missiles
    %Report
        SM_Vector = [SM_Fuel_Sub_WM(i),SM_Fuel_Sub_NM(i),SM_Fuel_Super_WM(i),SM_Fuel_Super_NM(i)];
        %Find upper bound point
        if placeholder1 == 0
            if all(SM_Vector <= 0.1)
                Upper_Limit(j) = d_display_WM(i);
                placeholder1 = 1;
            end
        end
        %Find lower bound point
        if placeholder2 == 0
            if any(SM_Vector <= -0.1)
                Lower_Limit(j) = d_display_NM(i);
                placeholder2 = 1;
            end
        end
        
    %Update d
        d_WM(i+1) = d_WM(i) + 1/Steps;
        d_NM(i+1) = d_NM(i) + 1/Steps;

        d_display_WM(i+1) = d_WM(i+1)*FL;
        d_display_NM(i+1) = d_NM(i+1)*FL;

end

%Update Fuel Weight
W_Fuel(j+1) = W_Fuel(j) + (1/Fuel_Step)*Wf;
FW_Display(j+1) = FW_Display(j) + (1/Fuel_Step)*100;
end

%Vector Adjustment
FW_Display(end) = [];

%plot
figure()
hold on
plot(FW_Display,Upper_Limit,'LineWidth',2)
plot(FW_Display,Lower_Limit,'LineWidth',2)
yline(0,'k')
xline(50,'--k')
xline(100*Wfu(end)/Wf,'--k')
fill([FW_Display, fliplr(FW_Display)],[Upper_Limit,fliplr(Lower_Limit)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
xlabel('Fuel Amount [%]')
ylabel('Distance from Empty CG [m]')
legend('Forewards Limit','Backwards Limit')

% Define the rectangle coordinates
x_rect = [0, 1, 1, 0]; % X coordinates
y_rect = [-10, -10, 4, 4]; % Y coordinates

% Create figure and fill the rectangle
hold on;
h = fill(x_rect, y_rect, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Match previous transparency
hold off;

% Remove from legend
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');


end