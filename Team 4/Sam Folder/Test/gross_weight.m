function [Wg_calc,EWF,Wf] = gross_weight(n_missile,n_eng,AR,T2W,WL,TSFC_cruise,TSFC_dash,L2D_cruise,L2D_dash,Ma_cruise,v_cruise,v_dash,a_dash,n_max,T2W_dash)

%% Standard Parameters
Ma = 1.6; %Design Mach Number
g = 9.81; %gravity

%% Payload

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

%% Start Numerical Solver
Wg(1) = 60000; %Initial Guess
error(1) = 20000; %Initial Guess
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

Wg_calc(end) = []; %Remove extra iteration

end