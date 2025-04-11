clear all; close all; clc;

%% Fixed Variables

rWp = 2.25*10^(-1); %Revenue per unit payload weight
rVp = 3.00*10^2; %Revenue per unit payload volume
c_e = 8.33*10^(-8); %Cost per unit energy
c_c = 5.00*10^(-1); %Cost (COC) per flight
cWg = 6.74*10^1; %Cost per unit gross weight
c_f = 5.56*10^(-5); %Cost (FOC) per unit flight-time
t_l = 3.60*10^6; %Total flight time in aircraft life

%% Selected Variables
EWF = 0.65; %empty weight fraction
FT = 3.5; %reference flight time [minutes]

%%
steps = 1000; %choose steps amount
Wg_min = 5; %minimum gross weight to test [N]
Wg_max = 50; %maximum gross weight to test [N]

%Convert reference time into proper units
Tf = FT*20*60; %seconds

Wg = linspace(Wg_min,Wg_max,steps);
for i = 1:steps
    %payload weight
        Wp(i) = Wg(i)*(1-EWF); %payload weight
        WpS(i) = (rWp*Wp(i))/(Tf); %payload weight score
    %gross weight
        WgS(i) = (cWg*Wg(i))/(t_l); %gross weight score
end

%% plot

figure()
hold on
plot(Wg,WpS)
plot(Wg,WgS)
legend('Payload Weight','Gross Weight','location','best')
xlabel('Gross Weight [N]')
ylabel('Term Score')