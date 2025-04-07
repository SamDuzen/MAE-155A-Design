clear all; close all; clc;

%% Parameter Definitions (Scoring Function)

rWp = 2.25*10^(-1); %Revenue per unit payload weight
rVp = 3.00*10^2; %Revenue per unit payload volume
c_e = 8.33*10^(-8); %Cost per unit energy
c_c = 5.00*10^(-1); %Cost (COC) per flight
cWg = 6.74*10^1; %Cost per unit gross weight
c_f = 5.56*10^(-5); %Cost (FOC) per unit flight-time
t_l = 3.60*10^6; %Total flight time in aircraft life

%% Scoring Parameter Initialization
Wp = 10; %Package 1 weight
Vp = 10; %Package 2 volume
Ef = 10; %Energy consumption
Wg = 10; %Gross weight
Tf = 10; %Flight time

%% Sensitivity Initialization
vary = -0.3:0.001:0.3;
J_base = ((rWp*Wp + rVp*Vp - c_e*Ef - c_c)/(Tf))-((cWg*Wg)/(t_l))-c_f;

%% Package 1 weight

for i = 1:length(vary)
    Wp_dummy = Wp + vary(i)*Wp; %percent change in Wp
    Sens_Wp(i) = ((rWp*Wp_dummy + rVp*Vp - c_e*Ef - c_c)/(Tf))-((cWg*Wg)/(t_l))-c_f;
    Sens_Wp_norm(i) = Sens_Wp(i)/J_base - 1;
end

%% Package 2 volume
for i = 1:length(vary)
    Vp_dummy = Vp + vary(i)*Vp; %percent change in Wp
    Sens_Vp(i) = ((rWp*Wp + rVp*Vp_dummy - c_e*Ef - c_c)/(Tf))-((cWg*Wg)/(t_l))-c_f;
    Sens_Vp_norm(i) = Sens_Vp(i)/J_base - 1;
end

%% Energy Consumption
for i = 1:length(vary)
    Ef_dummy = Ef + vary(i)*Ef; %percent change in Wp
    Sens_Ef(i) = ((rWp*Wp + rVp*Vp - c_e*Ef_dummy - c_c)/(Tf))-((cWg*Wg)/(t_l))-c_f;
    Sens_Ef_norm(i) = Sens_Ef(i)/J_base - 1;
end

%% Gross Weight
for i = 1:length(vary)
    Wg_dummy = Wg + vary(i)*Wg; %percent change in Wp
    Sens_Wg(i) = ((rWp*Wp + rVp*Vp - c_e*Ef - c_c)/(Tf))-((cWg*Wg_dummy)/(t_l))-c_f;
    Sens_Wg_norm(i) = Sens_Wg(i)/J_base - 1;
end

%% Flight Time
for i = 1:length(vary)
    Tf_dummy = Tf + vary(i)*Tf; %percent change in Wp
    Sens_Tf(i) = ((rWp*Wp + rVp*Vp - c_e*Ef - c_c)/(Tf_dummy))-((cWg*Wg)/(t_l))-c_f;
    Sens_Tf_norm(i) = Sens_Tf(i)/J_base - 1;
end

%% Plotting
figure()
hold on
plot(100*vary,100*Sens_Wp_norm)
plot(100*vary,100*Sens_Vp_norm)
plot(100*vary,100*Sens_Ef_norm)
plot(100*vary,100*Sens_Wg_norm)
plot(100*vary,100*Sens_Tf_norm)
legend('Package 1 Weight','Package 2 Volume','Energy Consumption','Gross Weight','Flight Time','Location','best')
xlabel('Percent Change in Scoring Parameter')
ylabel('Percent Change in Score')
grid on

