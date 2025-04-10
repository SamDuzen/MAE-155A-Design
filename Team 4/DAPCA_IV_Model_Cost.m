clear;close all;clc;
%% define input parameters
We = 2455.8; % empty weight [kg] ---> added up all of empty weights to and divided it by 9.81 m/s^2 to get kg 

FT = 6; %Flight time [hours]
Wf = 83536; %Fuel weight [N]
Fuel_Price = 711; %Price of JP-8 fuel per metric ton (1000kg) in USD (as of Feb 28, 2025)

Mach_max = 1.6;
speed_of_sound_mps = 343; % speed of sound at sea level (m/s) - adjust for altitude
V_mps = Mach_max * speed_of_sound_mps; % maximum velocity in m/s
V_kmh = V_mps * 3.6; % convert to km/h
V = V_kmh; % maximum velocity [km/h]


Q = 1000; % expected production quantity in five years 
FTA = 2; % number of flight-test aircraft

Neng = 2; % number of engines per aircraft

Tmax_lbf = 13500; % thrust per engine in lbf
Tmax_kN = Tmax_lbf * 4.448 / 1000; % total thrust in kN (converted from lbf)
Tmax = Tmax_kN; % engine maximum thrust [kN]

T_turbine_inlet = 1.6; % turbine inlet mach number
T_max_turbine_inlet = 1800; % turbine inlet temperature [K]
Cavionics = 1140400*Q; % avionics cost [dollars] (From AIAA document, avionics payload costs)


% choose the material and corresponding fudge factor
material = 'Aluminum'; % example, change as needed
switch material
    case 'Aluminum'
        fudge_factor = 1.0;
    case 'Graphite-epoxy'
        fudge_factor = 1.5; % taking the average of the range 1.1-1.8
    case 'Fiberglass'
        fudge_factor = 1.15; % taking the average of the range 1.1-1.2
    case 'Steel'
        fudge_factor = 1.75; % taking the average of the range 1.5-2.0
    case 'Titanium'
        fudge_factor = 1.45; % taking the average of the range 1.1-1.8
    otherwise
        error('Invalid material selected. Please choose from Aluminum, Graphite-epoxy, Fiberglass, Steel, or Titanium.');
end

%% DAPCA IV Cost Model

FF_General = 1.2; %General fudge factor for the DAPCA model for a more modern aircraft

% engineering hours (fps units)
He_fps = 4.86 * We^0.777 * V^0.894 * Q^0.163;

% engineering hours (mks units)
He_mks = 5.18 * We^0.777 * V^0.894 * Q^0.163;

% tooling hours (fps units)
Ht_fps = 5.99 * We^0.777 * V^0.696 * Q^0.263;

% tooling hours (mks units)
Ht_mks = 7.22 * We^0.777 * V^0.696 * Q^0.263;

% manufacturing hours (fps units)
Hm_fps = 7.37 * We^0.82 * V^0.484 * Q^0.641;

% manufacturing hours (mks units)
Hm_mks = 10.5 * We^0.82 * V^0.484 * Q^0.641;

% qc hours
HQ = 0.133 * Hm_mks; % using mks units for consistency

% development support cost (fps units)
Cd_fps = 91.3 * We^0.630 * V^1.3;

% development support cost (mks units)
Cd_mks = 67.4 * We^0.630 * V^1.3;

% flight test cost (fps units)
Cf_fps = 2498 * We^0.325 * V^0.822 * FTA^1.21;

% flight test cost (mks units)
Cf_mks = 1947 * We^0.325 * V^0.822 * FTA^1.21;

% manufacturing materials cost (fps units)
Cm_fps = 22.1 * We^0.921 * V^0.621 * Q^0.799;

% manufacturing materials cost (mks units)
Cm_mks = 31.2 * We^0.921 * V^0.621 * Q^0.799;

% engine production cost (fps units)
Ceng_fps = 3112/9.0667*Tmax + 243.25*max(1,T_turbine_inlet) + 0.9697*T_max_turbine_inlet -2228;

% engine production cost (mks units)
Ceng_mks = 3112/9.666*Tmax + 243.25*max(1,T_turbine_inlet) + 1.747*T_max_turbine_inlet -2228;
Ceng_mks = (245000000/59)*Q; %59 Ej200 engines were sold to the Spanish Air force for 245 Million dollars in Dec 24

%% fudge factor to relevant cost components
Cd_mks = Cd_mks * fudge_factor; % development Support Cost
Cm_mks = Cm_mks * fudge_factor; % manufacturing Materials Cost
%Ceng_mks = Ceng_mks * fudge_factor; % engine Production Cost

%% PRODUCTION COST 

% define production quantities for each year
production_quantities = [1:1000]; % FIX
learning_rate = 0.8; % 80% learning curve means 20% cost reduction per doubling of production)

% calculate production cost for each year using the learning curve
for i = 1:length(production_quantities)
  production_cost8(i) = Hm_mks * production_quantities(i)^(log10(learning_rate)/log10(2));
  production_cost9(i) = Hm_mks * production_quantities(i)^(log10(0.9)/log10(2));
  production_cost95(i) = Hm_mks * production_quantities(i)^(log10(0.95)/log10(2));
  production_cost7(i) = Hm_mks * production_quantities(i)^(log10(0.7)/log10(2));
  production_cost6(i) = Hm_mks * production_quantities(i)^(log10(0.6)/log10(2));
end



%% TOTAL 

% hourly wrap rates in 2024 (Converted from 2012)
Re = 157.12;
Rt = 161.22;
Rq = 147.56;
Rm = 133.9;

% total rdt&e + flyaway cost (using mks units)
RDTandE_flyaway = FF_General*(Cd_mks + Ht_mks*Rt + Cm_mks + Hm_mks*Rm + HQ*Rq + Cf_mks + 2*Ceng_mks + Cavionics + He_mks*Re);

%Flyaway cost 
C_flyaway = Ht_mks*Rt + Cm_mks + Hm_mks*Rm + Ceng_mks + Cavionics;


% display the results
disp(['Total RDT&E + Flyaway Cost: ', num2str(RDTandE_flyaway)]);

% Create a table using the calculated values
costTable = table(Cd_mks, Ht_mks, Cm_mks, Hm_mks, HQ, Cf_mks, Ceng_mks, Cavionics, RDTandE_flyaway, ...
    'VariableNames', { 'Development Support Cost', 'Tooling Hours Cost', 'Manufacturing Materials Cost', 'Manufacturing Hours Cost', ...
        'QC Hours Cost', 'Flight Test Cost', 'Engine Production Cost', 'Avionics Cost', 'Total RDT&E + Flyaway Cost' });

% Display the table
disp(costTable)


%% OPERATION AND SUPPORT COSTS (???) MAYBE GET RID OF 

% Define annual operating and support costs
annual_operating_costs = 1000000; % Example costs - adjust as needed
annual_support_costs = 5000000; % Example costs - adjust as needed

% Calculate total operating and support costs
total_operating_costs = sum(annual_operating_costs);
total_support_costs = sum(annual_support_costs);

% Fuel and Oil Costs
FHPY = 400; %Flight hours per year per aircraft
FPH = Wf/FT; %Fuel weight used per hour [N/hr]
FuelW_Yearly = FHPY*FPH; %Fuel weight used per year per aircraft [N/yr]
    Fuel_Total = (FuelW_Yearly/9.81)*(Fuel_Price/1000)*Q;

%Crew Salaries
CR = 1.1; %Crew Ratio
CN = 1; %Number of crew members per aircraft
CC = 2080*Re; %Cost per crew member
    Crew_Total = CN*CR*CC*Q; %Cost per year for crew to operate all produced aircrafts

% Maintenance
MMHPFH = 10; %Maintenance Man Hours Per Flight Hour
MMHPY = MMHPFH*FHPY; %Maintance Man Hours Per Year (per aircraft)
    Maintenance_Total = 2*MMHPY*Rm*Q; %Cost per year for maintance of all produced aircrafts

Operating_Total = Fuel_Total+Crew_Total+Maintenance_Total

%% TOTAL LIFE CYCLE

total_lifecycle_cost = RDTandE_flyaway + sum(production_cost8(1000)) + total_operating_costs + total_support_costs;

%% BREAK EVEN ANALYSIS 

% Define selling price per unit
selling_price = 25000000; % Example price - adjust as needed

% Calculate break-even quantity
break_even_quantity = (RDTandE_flyaway + total_operating_costs + total_support_costs) / abs((selling_price - Cm_mks/1000));

%% Display Results

costTable = table(RDTandE_flyaway, production_cost8(1000), total_operating_costs, total_support_costs, total_lifecycle_cost, ...
    'VariableNames', { 'RDT&E + Flyaway Cost', 'Production Costs', 'Operating Costs', 'Support Costs', 'Total Lifecycle Cost' });

disp(costTable)
disp(['Break-even Quantity: ', num2str(break_even_quantity)])



%% Visualization

figure;
hold on
plot(production_quantities, production_cost95);
plot(production_quantities, production_cost9);
plot(production_quantities, production_cost8);
plot(production_quantities, production_cost7);
plot(production_quantities, production_cost6);
xlabel('Production Quantity');
ylabel('Production Cost');
legend('95% Learning Rate', '90% Learning Rate', '80% Learning Rate','70% Learning Rate','60% Learning Rate')
grid on
title('Production Cost vs. Quantity');