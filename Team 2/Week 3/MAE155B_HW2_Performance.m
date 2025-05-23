%% Takeoff Analysis 

%Necessary values from Aerodynamics section
Clmax = 1.1066;
Wo = 22.2297;  %Gross weight [N]
Cdto = 0.01 ; %Drag force at 70% Vto and angle of attack = 0 [N]
%T =  ; %Thrust force at 70% Vto and angle of attack = 0 (From Prop section) [N]

% Approximated/Guessed values for now
WingLoading = 48*1.35; %Wing loading estimate from lecture, typical RC wing loading 1-3 lb/sqft, and conversion into N/m2 is multiplied by 48 [N/m^2]
S =0.2413; %Wing area by optimized gross weight and reference wing loading [m^2]
AR = 4.5; % Aspect ratio
e = 0.8; % Oswald efficiency

%Cd calculation
alpha0 = -2.669; %Zero lift angle of attack [deg]
dldalpha = 0.1182; % Rise over run approximation of 
Cl = 0.3; %Cl at 0 AoA from Sam's charts
Cdo = 0.01; % parasitic drag (zero lift) drag coefficient
K = 1/(pi*e*AR); 
Cd= 0.054; %Cdo + K*Cl^2;

%Constants
P = 101280; % Pressure in [Pa] 
Temp = 286 ; % Temperature in [K]
g = 9.81; % Gravitational constant [m/s^2]
Fc = 0.03; % Coefficient of rolling friction
RPM = 13000; %Propeller usage RPM 
P_diam = 9; %Propeller Diameter [in]
P_diam = P_diam*0.0254; %converting Prop Diam to meters [m]
rho = 1.23; % air density [kg/m^2]
n = RPM/60; %Prop revolutions per second
Ct = 0.11; %Coefficient of thrust at 70% Vto (38.3 mph)

% Calculations
Clto = Clmax*0.8; % Take off lift coefficient
Vto = sqrt(2*Wo/(Clto*rho*S)); 
L = 0.5*rho*(Clto)*S*(Vto*0.7)^2;
D = 0.5*rho*(Cd)*S*(Vto*0.7)^2;
T = rho* (n^2)* (P_diam^4)*Ct; %T = 0.998*9.81; %APC 9x4.5-E propeller Thrust [N]
amean = (g/Wo)*((T-D) - Fc*(Wo-L));
Sto = (Vto^2)/(2*amean); % Take off distance in [m]

fprintf('The Required takeoff distance is = %.3f m\n', Sto);

%% T vs V plot
K10 = readmatrix("APC 9x6E 13000 RPM data.xlsx");
J10000 = K10(:,2);
Ct10000 = K10(:,4);
Thrust = rho .* (n^2) .* (P_diam^4) .*Ct10000;
Velocity = P_diam .* n .* J10000;
V = linspace(3, 45, 421);
% % === Aerodynamic values from Seans workspace ===
% CL_cruise = 0.25535;
% CD_cruise = 0.0538775;
% K = 0.0011;
% CD0 = CD_cruise - K * CL_cruise^2;  % compute CD0 from drag polar
% 
% 
% % === DYNAMIC DRAG MODEL ===
% CL = (2 * Wo) ./ (rho * Velocity.^2 * S);
% CD = CD0 + K * CL.^2;

%Cl Calculation
CL = (2.*Wo)./(rho .* S .* Velocity.^2) ; % Coefficient of lift necessary to maintain lift at each V
K = 0.0011;
CL_min = -0.6174;
CD = 0.01 + K.*(CL-CL_min).^2 + + CL.^2/(pi*e*AR);

D = 0.5 * rho .* Velocity.^2 * S .* CD;
% O = length(Velocity);
% for k= 1:O
%     if D(k) > Thrust(k) && Velocity(K)<15
%         D(k) =      Thrust(k);
%     end
% end
% === FIND CRUISE POINT (T = D) ===
diff = Thrust - D;
idx_cross = find(diff(1:end-1) .* diff(2:end) < 0, 1,"last");
V1 = Velocity(idx_cross); V2 = Velocity(idx_cross+1);
T1 = Thrust(idx_cross); T2 = Thrust(idx_cross+1);
D1 = D(idx_cross); D2 = D(idx_cross+1);
V_cruise = V1 + (V2 - V1) * (T1 - D1) / ((T1 - D1) - (T2 - D2));
T_cruise = interp1(Velocity, Thrust, V_cruise);
nom_pct_thrust = 1.4389 / 11.5899 * 100;


figure;
plot(Velocity, Thrust, 'r-', 'LineWidth', 2); hold on;
plot(Velocity, D, 'b-', 'LineWidth', 2);hold on;
xline(24, 'k--', 'LineWidth', 2); hold on
plot(V_cruise, T_cruise, 'ko', 'MarkerFaceColor', 'g', 'DisplayName', 'Cruise Point');
xlabel('Velocity (m/s)');
ylabel('Force (N)');
ylim([0 25])
legend('Thrust', 'Drag','Cruise Speed', 'Maximum Speed point')
title('Thrust and Drag vs Velocity plot of APC 9x6E propeller')
grid on

pct_thrust = D./Thrust *100;
figure;
plot(Velocity,pct_thrust, 'k-', 'LineWidth', 2); hold on
plot(24,nom_pct_thrust,'ko', 'MarkerFaceColor', 'r')
xlim([10 max(Velocity)]);
ylim([0 100])
legend('Required thrust percentage','Nominal Cruise Speed point (12.41%)', Location='nw')
xlabel('Velocity (m/s)');
ylabel('Percentage of Thrust')
title('Percentage of maximum thrust needed for level flight');
grid on

%% Cl vs V plot 

% Key Velocities Calculations
Tmax = rho* (n^2)* (P_diam^4)*min(Ct10000);
H = length(Ct10000);
for k= 1:H
    if Ct10000(k) == min(Ct10000)
        Jmax = J10000(k);
        Vmax = P_diam*n*Jmax;
    end
end

Vstall = sqrt((2*Wo)/(rho*S*Clmax)); % Vstall from max Cl [m/s]
V = linspace(10, 45, 351);% Velocity range for plotting from 6 to 30 with 0.1 step size [m/s]


%Cl Calculation
CL = (2.*Wo)./(rho .* S .* V.^2) ; % Coefficient of lift necessary to maintain lift at each V

K = 0.0011;
CL_min = -0.6174;
CD = 0.01 + K.*(CL-CL_min).^2 + + CL.^2/(pi*e*AR);

% Plot
figure;
plot(V, CL, 'b-', 'LineWidth', 2); hold on;
xline(Vstall, 'r--','LineWidth', 2); hold on
xline(Vto, 'k--', 'LineWidth',2); hold on;
xline(V_cruise, 'g--', 'LineWidth',2); hold on;
xline(24, 'b--', 'LineWidth',2); 
xlabel('Velocity (m/s)');
ylabel('Coefficient of Lift');
legend('Cl', 'Stall Speed', 'Take off Speed', 'Maximum Speed','Nominal Cruise Speed');
title('Cl necessary to maintain steady level flight for varying speeds');
xlim([min(V) max(V)])
grid on;

