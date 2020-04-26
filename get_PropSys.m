function [rocket] = get_PropSys(rocket)

%% Propellant mass
% Chosen values
T = rocket.prop.T_avg;         %Thrust, [lbf], CHOSEN VALUE
Isp = rocket.prop.Isp;         %specific impulse, [s], CHOSEN VALUE
Itot = rocket.prop.I;          %total impulse, [lbf*s], GIVEN VALUE
OF = rocket.prop.OF;           %ox-fuel ratio, CHOSEN VALUE
P_c = rocket.prop.P_c;         %chamber pressure, [psf], CHOSEN VALUE
D_tank = rocket.geo.ox_t.D;    %tank diameter [ft]   CHOSEN VALUE

rocket.geo.fuel.D_i = 2/12;
D_i = rocket.geo.fuel.D_i;     %inner diam of fuel grain, [ft] CHOSEN VALUE

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

rho_htpb = 57.4;             %density cured htpb, [lb/ft3]
rho_paraffin = 56.2;         %density paraffin, [lb/ft3]

percent_htpb = .5;
percent_paraffin = .5;

%fuel density
rho_f = percent_htpb*rho_htpb + percent_paraffin*rho_paraffin; %[lb/ft3]

%oxidizer density
rho_ox = 48.21;            %density of liquid n2o [lb/ft3]

W_prop = Itot / Isp;           %propellant weight [lb]
rocket.weight.total.W_propellant = W_prop;

W_f = W_prop / (OF+1);         %fuel weight [lb]
W_ox = W_prop*OF / (OF+1);     %ox weight [lb]
rocket.weight.fuel.W = W_f;
rocket.weight.ox.W = W_ox;

%% Tank sizing

P_inj = 0.2 * P_c;              % [psf] pressure across injector, lower bound
P_feed = 50*144;                % [psf] pressure across feed sys
P_tank = P_c + P_inj + P_feed;  % [psf] ox tank pressure

V_ox = W_ox / rho_ox;       %ox volume [ft3]
L_tank = V_ox / (pi * D_tank^2 / 4);    %tank length [ft]
SA_tank = 2*(pi * D_tank^2 / 4) + pi*D_tank*L_tank;
t_tank = D_tank / 2 * P_tank / sig_working_tank;     %tank thickness, [ft]
rocket.geo.ox_t.L = L_tank + 4*t_tank; % approximate tank caps as 2t each
rocket.geo.ox_t.t = t_tank;
rocket.geo.ox.L = L_tank;
W_tank = SA_tank * t_tank * rho_tank; %tank mass, [lb]
rocket.weight.ox_t.W = W_tank;

%% Fuel grain + CC sizing
V_f = W_f / rho_f;          % [ft3]

D_cc = D_tank;
t_cc = D_cc/2 * P_c / sig_working_cc;   %assume cc is same material as tank
D_o = D_cc - 2*t_cc;               %outer diameter of fuel
rocket.geo.CC.D = D_cc;
rocket.geo.CC.t = t_cc;
rocket.geo.fuel.D_o = D_o;

L_f = V_f / (pi / 4 * (D_o^2 - D_i^2));     %fuel grain length, [m]
L_cc = L_f * 1.1;           %10% longer than fuel grain to allow mixing
rocket.geo.fuel.L = L_f;
rocket.geo.CC.L = L_cc;

SA_cc = pi*D_cc*L_cc;
W_cc = SA_cc * t_cc * rho_cc;       %cc mass [lb]
rocket.weight.CC.W = W_cc;

%nozzle
W_nozzle = 4*2.205;         %[lb]
d_nozzle_exit = 2.86/12;
A_e = pi * (d_nozzle_exit/2)^2;
rocket.prop.A_e = A_e;
%plumbing
%W_plumbing = 10*2.205;            %[lb]

rocket.weight.CC.W = rocket.weight.CC.W + W_nozzle;
%rocket.weight.misc.W = rocket.weight.misc.W + W_nozzle + W_plumbing;
