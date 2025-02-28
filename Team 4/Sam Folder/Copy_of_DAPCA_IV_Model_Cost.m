clear;close all;clc;
%% define input parameters
We = 2455.8; % empty weight [kg] ---> added up all of empty weights to and divided it by 9.81 m/s^2 to get kg 

YOO = 10; %Years of Operations

FT = 6.0284; %Flight time [hours]
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

% engineering hours (mks units)
He_mks = 5.18 * We^0.777 * V^0.894 * Q^0.163;

% tooling hours (mks units)
Ht_mks = 7.22 * We^0.777 * V^0.696 * Q^0.263;

% manufacturing hours (mks units)
Hm_mks = 10.5 * We^0.82 * V^0.484 * Q^0.641;

% qc hours
HQ = 0.133 * Hm_mks; % using mks units for consistency

% development support cost (mks units)
Cd_mks = 67.4 * We^0.630 * V^1.3;

% flight test cost (mks units)
Cf_mks = 1947 * We^0.325 * V^0.822 * FTA^1.21;

% manufacturing materials cost (mks units)
Cm_mks = 31.2 * We^0.921 * V^0.621 * Q^0.799;

% engine production cost (mks units)
Ceng_mks = (245000000/59)*Q*Neng; %59 Ej200 engines were sold to the Spanish Air force for 245 Million dollars in Dec 24

%% fudge factor to relevant cost components
Cd_mks = Cd_mks * fudge_factor; % development Support Cost
Cm_mks = Cm_mks * fudge_factor; % manufacturing Materials Cost
%Ceng_mks = Ceng_mks * fudge_factor; % engine Production Cost

%% TOTAL 

% hourly wrap rates in 2024 (Converted from 2012)
Re = 157.12;
Rt = 161.22;
Rq = 147.56;
Rm = 133.9;

% total rdt&e + flyaway cost (using mks units)
RDTandE_flyaway = FF_General*(Cd_mks + Ht_mks*Rt + Cm_mks + Hm_mks*Rm + HQ*Rq + Cf_mks + 2*Ceng_mks + Cavionics + He_mks*Re);

%% OPERATION AND SUPPORT COSTS (???) MAYBE GET RID OF 

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

Operating_Yearly = Fuel_Total+Crew_Total+Maintenance_Total; %Operating costs per year
Operating_Total = Operating_Yearly*YOO; %Lifetime operating costs

%% TOTAL LIFE CYCLE

%total_lifecycle_cost = RDTandE_flyaway + sum(production_cost8(1000)) + total_operating_costs + total_support_costs;
total_lifecycle_cost = RDTandE_flyaway + Operating_Total;

%% Display Results

summaryTable = table(RDTandE_flyaway, Operating_Yearly, Operating_Total, total_lifecycle_cost, ...
    'VariableNames', { 'RDT&E + Flyaway Cost', 'Yearly Operating Costs', 'Lifetime Operating Costs', 'Total Lifecycle Cost' });

operatingTable = table(Fuel_Total, Crew_Total, Maintenance_Total, Operating_Yearly, Operating_Total, ...
    'VariableNames', { 'Yearly Fuel Costs', 'Yearly Crew Costs', 'Yearly Maintenance Costs','Yearly Operating Costs','Lifetime Operating Costs' });

RDTable = table(He_mks*Re,Ht_mks*Rt,HQ*Rq,Hm_mks*Rm,Cd_mks,Cf_mks,Cm_mks,Ceng_mks, ...
    'VariableNames', { 'Engineering Costs [Total]','Tooling Costs [Total]','Quality Control Costs [Total]','Manufacturing Costs [Total]','Development Support Cost','Flight Tests Cost','Manufacturing Materials Cost','Engines Cost' });

disp(summaryTable)
disp(operatingTable)
disp(RDTable)

