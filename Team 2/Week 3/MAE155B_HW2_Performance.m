%% Takeoff Analysis 

%Necessary values from Aerodynamics section
Clmax = 1.5;
%S = ; %Wing area in  [m^2]
Wo = 13.6;  %Gross weight [N]
D = 5 ; %Drag force at 70% Vto and angle of attack = 0 [N]
%T =  ; %Thrust force at 70% Vto and angle of attack = 0 (From Prop section) [N]

% Approximated/Guessed values for now
WingLoading = 48*2.75; %Wing loading estimate from lecture, typical RC wing loading 1-3 lb/sqft, and conversion into N/m2 is multiplied by 48 [N/m^2]
S = Wo/WingLoading ; %Wing area by optimized gross weight and reference wing loading [m^2]
T = 0.998*9.81; %APC 9x4.5-E propeller Thrust [N]
L = Wo*0.8;

%Constants
P = 101280; % Pressure in [Pa] 
Temp = 286 ; % Temperature in [K]
rho = 1.23; % air density [kg/m^2]
g = 9.81; % Gravitational constant [m/s^2]
Fc = 0.03; % Coefficient of rolling friction

% Calculations
Clto = Clmax*0.8; % Take off lift coefficient
Vto = sqrt(2*Wo/(Clto*rho*S)); 
amean = (g/Wo)*((T-D) - Fc*(Wo-L));
Sto = (Vto^2)/(2*amean); % Take off distance in [m]

fprintf('The Required takeoff distance is = %.3f m\n', Sto);

%% T vs V plot
K10 = readmatrix("10000 RPM.xlsx");
J10000 = K10(:,2);
Ct10000 = K10(:,4);
RPM = 9860; %Propeller usage RPM 
P_diam = 9; %Propeller Diameter [in]
P_diam = P_diam*0.0254; %converting Prop Diam to meters [m]
rho = 1.23; % air density [kg/m^2]
n = RPM/60; %Prop revolutions per second
Thrust = rho .* (n^2) .* (P_diam^4) .*Ct10000;
Velocity = P_diam .* n .* J10000;
figure;
plot(Velocity, Thrust, 'r-', 'LineWidth', 2);
xlabel('Velocity (m/s)');
ylabel('Thrust (N)');
title('Thrust vs Velocity plot of APC 9x4.5 propeller')
grid on


%% Cl vs V plot 

Vstall = sqrt((2*Wo)/(rho*S*Clmax));
%Thrust vs V plot straight from propeller data 

% Velocity range (m/s)
V = linspace(10, 30, 201);  % speeds from 0 to 40 [m/s]

% Compute Lift and Drag
CL = (2.*Wo)./(rho .* S .* V.^2) ; % Coefficient of lift necessary to maintain
%steady level flight at specific Wo, S, for varying speeds
    %L = 0.5 * rho .* V.^2 * S .* CL;

% CD = CD0 + k .* CL.^2;
% D = 0.5 * rho .* V.^2 * S .* CD;


% Plot
figure;
plot(V, CL, 'b-', 'LineWidth', 2); hold on;

xlabel('Velocity (m/s)');
ylabel('Coefficient of Lift');
%legend('Lift', 'Drag');
title('Cl necessary to maintain steady level flight for varying speeds');
grid on;

