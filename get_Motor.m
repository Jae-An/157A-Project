function [rocket] = get_Motor(rocket)
%OUTPUT grain inner diam [IN],grain length [FT], mass flow rates [SLUGS/S],
%Jet Thrust [LBF], Isp [S],nozzle mass [slugs], Exit Area [ft^2], throat
%area [ft^2],expansion ratio
% ...propellant masses in SLUGS
% for an INPUT thrust [LBF], burn time {S], specific impulse [S]...
% chamber pressure [PSI], OF ratio, grain outer diam[IN]
% NOTE: returns zero inner grain diam if the input outer diam does not
% yield geometry with the required fuel mass flow rate 
g = 32.174;

%define parameters
T = rocket.prop.T_avg;         %Thrust, [lbf], CHOSEN VALUE
Isp = rocket.prop.Isp;         %specific impulse, [s], CHOSEN VALUE
Itot = rocket.prop.I;          %total impulse, [lbf*s], GIVEN VALUE
OF = rocket.prop.OF;           %ox-fuel ratio, CHOSEN VALUE
P_c = rocket.prop.P_c;         %chamber pressure, [psf], CHOSEN VALUE
D_oxt = rocket.geo.ox_t.D;     %tank diameter [ft]   CHOSEN VALUE
h_opt = rocket.prop.h_opt;     % altitude of optimization [ft]
% Material properties
FOS = 1.5;                     %factor of safety, CHOSEN VALUE
% Ox Tank
sig_yield_tank = rocket.geo.ox_t.material.yield;      %tensile strength, [psf]       
sig_working_tank = sig_yield_tank / FOS;              %working stress, [psf]
rho_tank = rocket.geo.ox_t.material.density;          %density, [lb/ft3]
% CC
sig_yield_cc = rocket.geo.CC.material.yield;          %tensile strength, [psf]       
sig_working_cc = sig_yield_cc / FOS;                  %working stress, [psf]
rho_cc = rocket.geo.CC.material.density;              %density, [lb/ft3]


%fuel density
rho_htpb = 1.785*g;             %density cured htpb, [lb/ft3]
rho_paraffin = 1.746*g;         %density paraffin, [lb/ft3]

percent_htpb = .5;
percent_paraffin = .5;

rho_f = percent_htpb*rho_htpb + percent_paraffin*rho_paraffin; %[lb/ft3]
rho_ox = 1.4984*g;                 %liquid n2o density [lb/ft3]

%calculate mass flow rate
Wdot_tot = T / Isp; %total mass flow rate, [lb/s]
Wdot_ox = OF / (OF + 1) * Wdot_tot;      %ox mass flowrate [lb/s]
Wdot_f_reqd = 1 / (OF + 1) * Wdot_tot;       %fuel mdot [lb/s]

%calculate weight
tb = Itot / T;
W_ox = Wdot_ox * tb;              %ox mass [sl]
W_f = Wdot_f_reqd * tb;         %fuel mass [sl]
theo_sliver_frac = 0.298;       %theoretical sliver fraction for circ. port
W_f = W_f * (1+theo_sliver_frac);   % account for leftover unburnt fuel

rocket.weight.fuel.W = W_f;
rocket.weight.ox.W = W_ox;
rocket.weight.total.W_propellant = W_f + W_ox;

%% Tank sizing
P_inj = 0.2 * P_c;              % [psf] pressure across injector, lower bound
P_feed = 50*144;                % [psf] pressure across feed sys
P_tank = P_c + P_inj + P_feed;  % [psf] ox tank pressure

V_ox = W_ox / rho_ox;       %ox volume [ft3]
L_tank = V_ox / (pi * D_oxt^2 / 4);    %tank length [ft]
SA_tank = 2*(pi * D_oxt^2 / 4) + pi*D_oxt*L_tank;
t_tank = D_oxt/2 * P_tank / sig_working_tank;     %tank thickness, [ft]
rocket.geo.ox_t.L = L_tank + 4*t_tank; % approximate tank caps as 2t each
rocket.geo.ox_t.t = t_tank;
rocket.geo.ox.L = L_tank;
W_tank = SA_tank * t_tank * rho_tank; %tank mass, [lb]
rocket.weight.ox_t.W = W_tank;

%% Size CC thickness + weight
V_f = W_f / rho_f;              %fuel volume [ft3]
D_cc = D_oxt;
t_cc = D_cc/2 * P_c / sig_working_cc;
D_o_fuel = D_cc - 2*t_cc;               %outer diameter of fuel
rocket.geo.CC.D = D_cc;
rocket.geo.CC.t = t_cc;
rocket.geo.fuel.D_o = D_o_fuel;

%% Fuel grain geometry
%initialize array of diameters
D_i_fuel = (0.1:0.1:D_o_fuel*12 - 0.1) / 12 ;  % inner diameter [ft]
t_init = 0.5 * (D_oxt - D_i_fuel);     % initial thickness of grain, [ft]

% calc grain geometries
L_grain = V_f ./ (pi/4 * (D_o_fuel^2 - D_i_fuel.^2)); %grain length, [ft]
A_cs = pi ./ 4 * D_i_fuel.^2;         % grain port area   [ft2]
A_b = pi * D_i_fuel .* L_grain;       %grain burn area  [ft2]

LD = L_grain ./ D_oxt;            %length-diam ratio of grain
%stanford hybrid paper calls for L/D >= 4, for proper propellant mixing (?)

%% regression rate
a = 0.1146;                     %from paraffin/htpb/n2o research paper
n = 0.5036;                     % for G_o in units of [lbm/in2/s]

G_o = (Wdot_ox/g) ./ A_cs;           %oxidizer mass velocity [sl/ft2/s]
G_o = G_o * g / 144;         %ox mass vel. [lbm/in2/s]

rdot = a * G_o.^n;              %regression rate [in/s]
rdot = rdot / 12;               %reg rate [ft/s]
Wdot_f = rho_f * A_b .* rdot;   % fuel mass flow rate [lb/s]

grain_d_i = interp1(Wdot_f,D_i_fuel, Wdot_f_reqd,'linear',0); % [ft]
grain_length = V_f ./ (pi/4 * (D_oxt^2 - grain_d_i.^2)); % [ft]
rocket.geo.fuel.D_i = grain_d_i;
rocket.geo.fuel.L = grain_length;
rocket.geo.CC.L = grain_length*1.1;

Ac_cc = pi/4 * (D_cc^2 - (D_cc-2*t_cc)^2);
W_cc = Ac_cc * rocket.geo.CC.L * rho_cc;       %cc mass [lb]
rocket.weight.CC.W = W_cc;

%% NOZZLE CALCULATIONS
%constants
R = 8.314; %universal gas constant (J/mol*K)
[~,P_e,~] = get_Atmosphere(0,h_opt);

%design choices
T_c = 1500; %(K)                                      %CHOICE                   
mdot = Wdot_tot / g; %[slug/s]                          %FROM Above

%combustion chamber parameters
R_specific = R/0.030;   %[J/kg K]                     % from combustion_products.m
R_specific = R_specific*0.737562/0.0685218; %[ft*lb/slug K]
k = 1.25111;                                 %from proprep3 program

%Conical Nozzle assumption
alpha = 15; %(degrees), Divergent cone half angle      %CHOSEN
alpha = alpha*pi/180; %Radians
cone_thickness = 1/100 *3.281; %ft %thickness of cone walls    %Guess
rho_nozzle = 1.8*100^3 /1000; %graphite density (kg/m^3)     %choice of material
rho_nozzle = rho_nozzle/515; %slug/ft^3
lambda = (1 + cos(alpha))/2;

%throat calculations
P_t = P_c*(1+ (k-1)/2)^(-k/(k+1)); %throat pressure
T_t = T_c/(1 + (k-1)/2); %throat temperature
A_t = (mdot/P_t)*sqrt(R_specific*T_t); %area of throat
R_t = sqrt(A_t/pi); %radius of throat

%Nozzle calculations (more)
Ma_e = ((2/(k-1))*((P_c/P_e)^((k-1)/k) - 1))^0.5; %exit mach number
A_e = (A_t/Ma_e)*((2 + Ma_e^2 *(k-1))/(k+1))^((k+1)/(2*(k-1))); %exit area (ft^2)
R_e = sqrt(A_e/pi); %exit radius (ft)
L_n = (R_e - R_t)/tan(alpha); %cone length (ft)
V_inside_cone = (1/3)*pi*L_n*(R_t^2 + R_t*R_e + R_e^2);  %(ft^3)
R_e_out = R_e + cone_thickness;
R_t_out = R_t + cone_thickness;
V_outside_cone = (1/3)*pi*L_n*(R_t_out^2 + R_t_out*R_e_out + R_e_out^2);
V_cone = V_outside_cone - V_inside_cone; %volume of nozzle
nozzle_mass = rho_nozzle*V_cone; %mass of nozzle (slug)
epsilon = A_e/A_t; %expansion ratio

rocket.prop.A_e = A_e; % doesn't include nozzle wall thickness
rocket.geo.misc.L(7) = L_n;
rocket.geo.misc.x(7) = rocket.geo.total.L - L_n;
rocket.weight.misc.W(7) = nozzle_mass * g;
rocket.weight.misc.CG(7) = rocket.geo.misc.x(7) + L_n/2;

%Exit velocity
u_e = sqrt(2*k*R_specific*T_c/(k-1) *(1 - P_e/P_c)^((k-1)/k)); %exit velocity (ft/s)
%textbook value of u_e for HTPB/n20 at 1000 psi P_c = 1625 m/s

%thrust
Thrust_jet = lambda*mdot*u_e;
rocket.prop.T_avg = Thrust_jet;

%max launch weight
%Max_weight = Thrust_launch/(V_rail^2 /(2*L_rail*g) + 1); %in lbf
Isp_new = u_e/g; %specific impulse
rocket.prop.Isp = Isp_new;

end