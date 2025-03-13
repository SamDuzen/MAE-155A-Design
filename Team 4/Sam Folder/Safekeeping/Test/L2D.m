function [L2D_cruise,L2D_dash,L2D_combat,v_cruise,v_dash,a_dash] = L2D(Ma_cruise,Cdo_c,Cdo_d,WL_Cruise,WL_dash,AR,e,e_combat,n_max)

%% Global Parameters
Ma = 1.6; %Design Mach Number

%% Improved L/D Estimate (For cruise/loiter)

%Dynamic Pressure @ Cruise - 9500 m (31000 ft)
den_95 = 0.43966; %[kg/m^3]
a_95 = 301.7; %[m/s]
v_cruise = Ma_cruise*a_95; %cruise velocity
q_cruise = 0.5*den_95*v_cruise^2; %dynamic pressure [N/m^2]

%L/D Estimation
L2D_cruise = 1/(((q_cruise*Cdo_c)/(WL_Cruise))+WL_Cruise*(1/(q_cruise*pi*AR*e)));

%% L/D Estimate (For dash)
    %no known better method yet so for now...
a_dash = 297.4; %speed of sound at 10.5km (about 35000 ft) [m/s]
v_dash = Ma*a_dash;
den_dash = 0.388857; %[kg/m^3]
q_dash = 0.5*den_dash*v_dash^2;

L2D_dash = 1/(((q_dash*Cdo_d)/(WL_dash))+WL_dash*(1/(q_dash*pi*AR*e)));

%% L/D Estimate (for combat)
L2D_combat = 1/(((q_dash*Cdo_d)/(n_max*WL_dash))+WL_dash*n_max*(1/(q_dash*pi*AR*e_combat)));

end