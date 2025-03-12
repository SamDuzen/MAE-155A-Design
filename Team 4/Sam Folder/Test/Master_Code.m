clear all; close all; clc;

%% Input Parameters

%Currently all input parameters are in the "Parameters.tsv" file. Make
%adjustments there.

%% TSV File Data Input - Input Parameters

t = readtable("Parameters.tsv", "FileType","text",'Delimiter', '\t');
    n_missile = t.(2)(1);
    n_eng = t.(2)(2);
    Ma_cruise = t.(2)(3);
    e = t.(2)(4);
    e_combat = t.(2)(5);
    n_max = t.(2)(6);
    WL = t.(2)(7);
    WL_Cruise = t.(2)(8);
    WL_dash = t.(2)(9);
    T2W = t.(2)(10);
    T2W_dash = t.(2)(11);
    TSFC = t.(2)(12);
    TSFC_cruise = t.(2)(13);
    TSFC_dash = t.(2)(14);
    Cdo_c = t.(2)(15);
    Cdo_d = t.(2)(16);

%% Initialize
Wg = 122917.9; %gross weight [N]


%% Run Sub Functions

%Run climb angle calculations
Ma_test = 0.9;
[CA] = ClimbAngle(Ma_test);

%Run wing sizing code to get wing parameters
[Wing_Table,S,b,C_Root,C_Tip,MAC,Y_MAC,Sw_LE,Sw_QC,TaR,ThR,AR,WTA,WIA,DA,S_VT,S_HT,b_HT,C_Root_HT,C_Tip_HT,AR_HT,TaR_HT,MAC_HT,Y_MAC_HT,MAC_VT,Y_MAC_VT,b_VT,C_Root_VT,C_Tip_VT,AR_VT,TaR_VT] = wing_size(Wg);

%Run L2D function
[L2D_cruise,L2D_dash,L2D_combat,v_cruise,v_dash,a_dash] = L2D(Ma_cruise,Cdo_c,Cdo_d,WL_Cruise,WL_dash,AR,e,e_combat,n_max);

%Run gross weight estimation function
[Wg_calc,EWF,Wf,Flight_Time] = gross_weight(n_missile,n_eng,AR,T2W,WL,TSFC_cruise,TSFC_dash,L2D_cruise,L2D_dash,Ma_cruise,v_cruise,v_dash,a_dash,n_max,T2W_dash);
    Wg = Wg_calc(end); %Save last vector entry as converged gross weight

%Run weight balance function
[CG_WM,CG_NM,X_CG_WM,X_CG_NM,NP_Sub,NP_Super,SM_Sub_WM,SM_Sub_NM,SM_Super_WM,SM_Super_NM,d_display] = Weight_Balance(Wg,Ma_cruise,Wf);

%Find Fuel Volume
Fuel_Volume = (Wf/9.81)/119.826;

