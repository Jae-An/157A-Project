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
tb = Itot/T;
OF = rocket.prop.OF;           %ox-fuel ratio, CHOSEN VALUE
P_c = rocket.prop.P_c / 144;         %chamber pressure, [psi], CHOSEN VALUE
D_oxt = 0.5 * 12;%rocket.geo.ox_t.D;     % [in]tank diameter [in]   CHOSEN VALUE
h_opt = rocket.prop.h_opt;     % altitude of optimization [ft]
% Material properties
FOS = 1.5;                     %factor of safety, CHOSEN VALUE
% Ox Tank
sig_yield_tank = rocket.geo.ox_t.material.yield;      %tensile strength, [psf]       
sig_working_tank = sig_yield_tank / FOS;              %working stress, [psf]
rho_tank = rocket.geo.ox_t.material.density / g;          %density, [slug/ft3]
% CC
sig_yield_cc = rocket.geo.CC.material.yield;          %tensile strength, [psf]       
sig_working_cc = sig_yield_cc / FOS;                  %working stress, [psf]
rho_cc = rocket.geo.CC.material.density / g;              %density, [slug/ft3]


%fuel density
rho_htpb = 1.785;             %density cured htpb, [sl/ft3]
rho_paraffin = 1.746;         %density paraffin, [sl/ft3]

percent_htpb = .5;
percent_paraffin = .5;

rho_f = percent_htpb*rho_htpb + percent_paraffin*rho_paraffin; %[lb/ft3]
rho_ox = 1.4984;                 %liquid n2o density [sl/ft3]

%calculate mass flow rate
mdot_tot = T / Isp /g; %total mass flow rate, [sl/s]
mdot_ox = OF / (OF + 1) * mdot_tot;      %ox mass flowrate [sl/s]
mdot_f_reqd = 1 / (OF + 1) * mdot_tot;       %fuel mdot [sl/s]

%calculate weight
m_ox = mdot_ox * tb;              %ox mass [sl]
m_f = mdot_f_reqd * tb;         %fuel mass [sl]
theo_sliver_frac = 0.298;       %theoretical sliver fraction for circ. port
m_f = m_f * (1+theo_sliver_frac);   % account for leftover unburnt fuel

rocket.weight.fuel.W = m_f * g;
rocket.weight.ox.W = m_ox * g;
rocket.weight.total.W_propellant = (m_f + m_ox)*g;

%% Tank sizing
P_inj = 0.2 * P_c;              % [psi] pressure across injector, lower bound
P_feed = 50;                % [psi] pressure across feed sys
P_tank = (P_c + P_inj + P_feed) * 144;  % [psf] ox tank pressure
t_tank = (D_oxt/12)/2 * P_tank / sig_working_tank;     % [ft]tank thickness, [ft]

V_ox = m_ox / rho_ox;       %ox volume [ft3]
L_tank = (V_ox-((pi/6)*((D_oxt/12)-2*t_tank)^3)) / (pi * ((D_oxt/12)-2*t_tank)^2 / 4); %cylindrical tank length [ft]
V_tank = (pi/6)*((D_oxt/12)^3 - ((D_oxt/12)-2*t_tank)^3) + (pi/4)*L_tank*((D_oxt/12)^2 - ((D_oxt/12)-2*t_tank)^2);
W_tank = V_tank * rho_tank * g; %tank mass, [lb]

rocket.geo.ox_t.L = L_tank + D_oxt/12; % assume hemi caps
rocket.geo.ox_t.t = t_tank;
rocket.geo.ox.L = L_tank;
rocket.weight.ox_t.W = W_tank;

%% Size CC thickness + weight
V_f = m_f / rho_f;              %fuel volume [ft3]
D_cc = D_oxt/12;
t_cc = D_cc/2 * (P_c*144) / sig_working_cc;
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

G_o = (mdot_ox) ./ A_cs;           %oxidizer mass velocity [sl/ft2/s]
G_o = G_o * g / 144;         %ox mass vel. [lbm/in2/s]

rdot = a * G_o.^n;              %regression rate [in/s]
rdot = rdot / 12;               %reg rate [ft/s]
mdot_f = rho_f * A_b .* rdot;   % fuel mass flow rate [lb/s]

grain_d_i = interp1(mdot_f,D_i_fuel, mdot_f_reqd,'linear',0); % [ft]
grain_length = V_f ./ (pi/4 * (D_o_fuel^2 - grain_d_i.^2)); % [ft]
rocket.geo.fuel.D_i = grain_d_i;
rocket.geo.fuel.L = grain_length;
rocket.geo.CC.L = grain_length*1.1;

Ac_cc = pi/4 * (D_cc^2 - (D_cc-2*t_cc)^2);
W_cc = Ac_cc * rocket.geo.CC.L * rho_cc * g;       %cc mass [lb]
rocket.weight.CC.W = W_cc;

%% NOZZLE CALCULATIONS
%constants
R = 8.314; %universal gas constant (J/mol*K)
[~,P_e,~] = get_Atmosphere(0,h_opt);
P_e = P_e / 144; % [psi]

%design choices
mdot = mdot_tot; %[slug/s]                          %FROM Above

%combustion chamber parameters
R_specific = R/0.025427;   %[J/kg K]                     % from combustion_products.m
R_specific = R_specific*0.737562/0.0685218; %[ft*lb/slug K]
k = 1.25111;                                 %from proprep3 program

%combustion efficiency and chamber temp
cstar_theo = 5232.68;   % theoretical characteristic velocity, [ft/s] from proprep3
comb_eff = .80;         % combustion efficiency
cstar = cstar_theo * comb_eff;  %char. vel., [ft/s]
T_c = (cstar/(2/(k+1))^(-(k+1)/(2*(k-1))))^2 /(R_specific/k);

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
A_e = (A_t/Ma_e)*((2 + Ma_e^2 *(k-1))/(k+1))^((k+1)/(2*(k-1))); %exit area (in^2)
R_e = sqrt(A_e/pi); %exit radius (in)
L_n = (R_e - R_t)/tan(alpha); %cone length (in)
V_inside_cone = (1/3)*pi*L_n*(R_t^2 + R_t*R_e + R_e^2);  %(in^3)
R_e_out = R_e + cone_thickness;
R_t_out = R_t + cone_thickness;
V_outside_cone = (1/3)*pi*L_n*(R_t_out^2 + R_t_out*R_e_out + R_e_out^2);
V_cone = V_outside_cone - V_inside_cone; %volume of nozzle
nozzle_mass = rho_nozzle * (V_cone / 12^3); %mass of nozzle (slug)
epsilon = A_e/A_t; %expansion ratio

rocket.prop.A_e = A_e / 144; % doesn't include nozzle wall thickness
rocket.geo.misc.L(8) = L_n / 12;
rocket.geo.misc.x(8) = rocket.geo.total.L - (L_n/12);
rocket.weight.misc.W(8) = nozzle_mass * g;
rocket.weight.misc.CG(8) = rocket.geo.misc.x(7) + (L_n/12)/2;

%Exit velocity
u_e = sqrt(2*k*R_specific*T_c/(k-1) *(1 - (P_e/P_c)^((k-1)/k))); %exit velocity (ft/s)
%textbook value of u_e for HTPB/n20 at 1000 psi P_c = 1625 m/s

%thrust
Thrust_jet = lambda*mdot*u_e;
rocket.prop.T_avg = Thrust_jet;

%max launch weight
%Max_weight = Thrust_launch/(V_rail^2 /(2*L_rail*g) + 1); %in lbf
Isp_new = u_e/g; %specific impulse
rocket.prop.Isp = Isp_new;

end