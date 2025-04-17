%% Takeoff Analysis 

%Necessary values from Aerodynamics section
Clmax = 1.47;
Wo = 13.6;  %Gross weight [N]
Cdto = 0.01 ; %Drag force at 70% Vto and angle of attack = 0 [N]
%T =  ; %Thrust force at 70% Vto and angle of attack = 0 (From Prop section) [N]

% Approximated/Guessed values for now
WingLoading = 48*1.35; %Wing loading estimate from lecture, typical RC wing loading 1-3 lb/sqft, and conversion into N/m2 is multiplied by 48 [N/m^2]
S = Wo/WingLoading ; %Wing area by optimized gross weight and reference wing loading [m^2]
AR = 5; % Aspect ratio
e = 0.8; % Oswald efficiency

%Cd calculation
alpha0 = -6; %Zero lift angle of attack [deg]
dldalpha = 1.2/12; % Rise over run approximation of 
Cl = 0.6; %Cl at 0 AoA from Sam's charts
Cdo = 0.01; % parasitic drag (zero lift) drag coefficient
K = 1/(pi*e*AR); 
Cd= Cdo + K*Cl^2;

%Constants
P = 101280; % Pressure in [Pa] 
Temp = 286 ; % Temperature in [K]
g = 9.81; % Gravitational constant [m/s^2]
Fc = 0.03; % Coefficient of rolling friction
RPM = 9860; %Propeller usage RPM 
P_diam = 9; %Propeller Diameter [in]
P_diam = P_diam*0.0254; %converting Prop Diam to meters [m]
rho = 1.23; % air density [kg/m^2]
n = RPM/60; %Prop revolutions per second
Ct = 0.042; %Coefficient of thrust at 70% Vto (38.3 mph)

% Calculations
Clto = Clmax*0.8; % Take off lift coefficient
Vto = sqrt(2*Wo/(Clto*rho*S)); 
L = 0.5*rho*(Clto)*S*(Vto*0.7)^2;
D = 0.5*rho*(Cd)*S*(Vto*0.7)^2;
T = rho* (n^2)* (P_diam^4)*Ct ;%T = 0.998*9.81; %APC 9x4.5-E propeller Thrust [N]
amean = (g/Wo)*((T-D) - Fc*(Wo-L));
Sto = (Vto^2)/(2*amean); % Take off distance in [m]

fprintf('The Required takeoff distance is = %.3f m\n', Sto);

%% T vs V plot
K10 = readmatrix("10000 RPM.xlsx");
J10000 = K10(:,2);
Ct10000 = K10(:,4);
Thrust = rho .* (n^2) .* (P_diam^4) .*Ct10000;
Velocity = P_diam .* n .* J10000;
figure;
plot(Velocity, Thrust, 'r-', 'LineWidth', 2);
xlabel('Velocity (m/s)');
ylabel('Thrust (N)');
title('Thrust vs Velocity plot of APC 9x4.5 propeller')
grid on


%% Cl vs V plot 

% Key Velocities Calculations
Tmax = rho* (n^2)* (P_diam^4)*min(Ct10000);
H = length(Ct10000);
for k= 1:H
    if Ct10000(k) == 0
        Jmax = J10000(k);
        Vmax = P_diam*n*Jmax;
    end
end
Vstall = sqrt((2*Wo)/(rho*S*Clmax)); % Vstall from max Cl [m/s]
V = linspace(6, 30, 241);% Velocity range for plotting from 6 to 30 with 0.1 step size [m/s]


%Cl Calculation
CL = (2.*Wo)./(rho .* S .* V.^2) ; % Coefficient of lift necessary to maintain lift at each V


% Plot
figure;
plot(V, CL, 'b-', 'LineWidth', 2); hold on;
xline(Vstall, 'r--','LineWidth', 2); hold on
xline(Vto, 'k--', 'LineWidth',2);
xline(Vmax, 'g--', 'LineWidth',2);
xlabel('Velocity (m/s)');
ylabel('Coefficient of Lift');
legend('Cl', 'Stall Speed', 'Take off Speed', 'Maximum Speed');
title('Cl necessary to maintain steady level flight for varying speeds');
xlim([min(V) max(V)])
grid minor;

