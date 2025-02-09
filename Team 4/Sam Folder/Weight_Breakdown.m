clear all; close all; clc;

%%

% THESE ARE ALL IN BU UNITS SO WHEN IMPLEMENTING THE NUMBER, CONVERT TO BU,
% ILL MAKE A CODE LATER OR WE CAN JUST LOOK ONLINE BUT FOR NOW REMEMEBER
% ITS ALL IN BU UNITS 


%% Adjustable Parameters

Kdw = 1; % delta wing factor (1 for non-delta wings)
Kvs = 1; % variable sweep factor (1 for fixed wings)
Wdg = 100000; % flight design gross weight (lb)
Nz = 3.5; % ultimate load factor (limit load factor x 1.5)
Sw = 300; % trapezoidal wing area (ft^2)
A = 9; % wing aspect ratio
tc_root = 0.12; % thickness-to-chord ratio at the root
lambda = 1; % taper ratio of the wing (root chord/tip chord)
Lambda = 25; % wing sweep angle at 25% mac (degrees)
Scsw = 5; % wing-mounted control surface area (ft^2)

% fuselage and structure parameters
Fw = 10; % fuselage width at horizontal tail intersection (ft)
Bh = 4; % fuselage structural depth (ft)
Sht = 50; % horizontal tail area (ft^2)
Krht = 1.1; % rolling horizontal tail factor
Ht = 2; % horizontal tail height above fuselage (ft)
Hv = 5; % vertical tail height above fuselage (ft)
Svt = 30; % vertical tail area (ft^2)
M = 0.85; % mach number (design maximum)
Lt = 15; % tail length from wing quarter-mac to tail quarter-mac (ft)
Sr = 3; % rudder area (ft^2)
Avt = 2; % vertical tail aspect ratio
Lambda_vt = 20; % vertical tail sweep angle (degrees)

% miscellaneous structural factors
Kdwf = 1.05; % delta-wing fuselage factor
L = 60; % fuselage structural length (ft)
D = 15; % engine diameter (ft)
W = 80000; % total fuselage structural weight (lb)
Kcb = 1.02; % cross-beam factor for gear
Ktpg = 1.03; % tripod gear factor
Wl = 7000; % landing design gross weight (lb)
Nl = 2; % ultimate landing load factor (ngear x 1.5)

% landing gear parameters
Lm = 4; % extended length of main landing gear (in)
Ln = 5; % extended nose gear length (in)
Nnw = 2; % number of nose wheels
Nen = 2; % number of engines
T = 20000; % total engine thrust (lb)

% engine parameters
Sfw = 20; % firewall surface area (ft^2)
Wen = 5000; % engine weight, each (lb)
Kvg = 1.1; % variable geometry factor for air induction
Ld = 2; % duct length (ft)
Kd = 0.9; % duct constant (dimensionless)
Ls = 0.5; % single duct length (ft)
De = 1.2; % engine diameter (ft)
Ltp = 4; % length of tailpipe (ft)
Lsh = 3; % length of engine section shell (ft)
Lec = 2.5; % routing distance from engine front to cockpit (ft)
Te = 18000; % thrust per engine (lb)

% fuel parameters
Vt = 1000; % total fuel volume (gal)
Vi = 200; % integral tank volume (gal)
Vp = 300; % self-sealing "protected" tanks volume (gal)
Nt = 2; % number of fuel tanks
SFC = 0.5; % specific fuel consumption at maximum thrust (lb/hr/lb)

% crew and control parameters
Nc = 1; % number of crew
Scs = 200; % total area of control surfaces (ft^2)
Ns = 2; % number of flight control systems
Nci = 1; % number of crew equivalents (1.0 for single pilot)
Kvsh = 1.05; % variable sweep wing hydraulic system factor
Nu = 4; % number of hydraulic utility functions
Kmc = 1.1; % mission completion factor
Rkva = 300; % system electrical rating (kv·a)
La = 20; % electrical routing distance (ft)
Ngen = 2; % number of generators
Wuav = 25000; % uninstalled avionics weight (lb)


%% Calculations
W_wing = 0.0103 * Kdw * Kvs * (Wdg * Nz)^0.5 * Sw^0.622 * A^0.785 * (tc_root)^0.785 * (1 + lambda)^0.05 * (cosd(Lambda))^(-1) * Scsw^0.04;

W_horizontal_tail = 3.316 * (1 + Fw / Bh)^(-2) * (Wdg * Nz / 1000)^0.26 * Sht^0.806;

W_vertical_tail = 0.452 * Krht * (1 + Ht / Hv)^0.5 * (Wdg * Nz)^0.488 * Svt^0.718 * M^0.341 * Lt^(-1) * (1 + Sr / Svt)^0.348 * Avt^0.223 * ...
                   (1 + lambda)^0.25 * (cosd(Lambda_vt))^(-0.323);

W_fuselage = 0.499 * Kdwf * Wdg^0.35 * Nz^0.25 * L^0.5 * D^0.849 * W^0.685;

W_main_landing_gear = Kcb * Ktpg * (Wl * Nl)^0.25 * Lm^0.973;

W_nose_landing_gear = (Wl * Nl)^0.29 * Ln^0.5 * Nnw^0.525;

W_engine_mounts = 0.013 * Nen^0.795 * T^0.579 * Nz;

W_firewall = 1.13 * Sfw;

W_engine_section = 0.01 * Wen^0.717 * Nen * Nz;

% If Ls/Ld < 0.25, set ratio to 0.25
ratio_Ls_Ld = max(Ls/Ld, 0.25);
W_air_induction_system = 13.29 * Kvg * Ld^0.643 * Kd^0.182 * Nen^1.498 * (ratio_Ls_Ld)^(-0.373) * De;

W_tailpipe = 3.5 * De * Ltp * Nen;

W_engine_cooling = 4.55 * De * Lsh * Nen;

W_oil_cooling = 37.82 * Nen^1.023;

W_engine_controls = 10.5 * Nen^1.008 * Lec^0.222;

W_starter_pneumatic = 0.025 * Te^0.760 * Nen^0.72;

W_fuel_system_tanks = 7.45 * Vt^0.47 * (1 + Vi/Vt)^(-0.095) * (1 + Vp/Vt) * Nt^0.066 * Nen^0.052 * (T * SFC / 1000)^0.249;

W_flight_controls = 36.28 * M^0.003 * Scs^0.489 * Ns^0.484 * Nc^0.127;

W_instruments = 8.0 + 36.37 * Nen^0.676 * Nt^0.237 + 26.4 * (1 + Nci)^1.356;

W_hydraulics = 37.23 * Kvsh * Nu^0.664;

W_electrical = 172.2 * Kmc * Rkva^0.152 * Nc^0.10 * La^0.10 * Ngen^0.091;

W_avionics = 2.117 * Wuav^0.933;

W_furnishings = 217.6 * Nc; % includes seats

W_air_conditioning_antiice = 201.6 * ((Wuav + 200 * Nc) / 1000)^0.735;

W_handling_gear = 3.2e-4 * Wdg;

%% Creating Table
component_names = {'Wing', 'Horizontal Tail', 'Vertical Tail', 'Fuselage', 'Main Landing Gear', ...
    'Nose Landing Gear', 'Engine Mounts', 'Firewall', 'Engine Section', 'Air Induction System', ...
    'Tailpipe', 'Engine Cooling', 'Oil Cooling', 'Engine Controls', 'Starter Pneumatic', ...
    'Fuel System Tanks', 'Flight Controls', 'Instruments', 'Hydraulics', 'Electrical System', ...
    'Avionics', 'Furnishings', 'Air Conditioning & Anti-Ice', 'Handling Gear'};

weights = [W_wing, W_horizontal_tail, W_vertical_tail, W_fuselage, W_main_landing_gear, ...
    W_nose_landing_gear, W_engine_mounts, W_firewall, W_engine_section, W_air_induction_system, ...
    W_tailpipe, W_engine_cooling, W_oil_cooling, W_engine_controls, W_starter_pneumatic, ...
    W_fuel_system_tanks, W_flight_controls, W_instruments, W_hydraulics, W_electrical, ...
    W_avionics, W_furnishings, W_air_conditioning_antiice, W_handling_gear];

% display as  table
T = table(component_names', weights', 'VariableNames', {'Component', 'Weight_lb'});
disp(T);
