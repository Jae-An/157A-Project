function [t_tank, L_tank,  m_tank] = get_tank(p_c,m_o, d_tank_o)
%sizes tank
% input chamber pressure in PSI, ox mass in SLUGS, d_tank in INCHES
% output tank thickness in INCHES, tank length in FEET, tank mass in POUNDS
%NOTE: does not take into account ullage space


%inputs
d_tank_o = d_tank_o /12;               %tank diam [ft]

p_c = p_c * 144;                %chamber pressure, [psf]
p_inj = 0.2 * p_c;              %pressure across injector, lower bound
p_feed = 50 * 144;      %pressure across feed sys [psf]
p_tank = p_c + p_inj + p_feed;  %ox tank pressure [psf]

%Al 7050-T7651 properties
sig_yield = 1.0234 * 10^7;          %tensile strength, [psf]
FOS = 1.5;                  %factor of safety, CHOSEN VALUE
sig_working = sig_yield / FOS;  %working stress, [psf]
rho_Al = 5.4911;              %density, [sl/ft3]

%calc tank dimensions
t_tank = d_tank_o / 2 * p_tank / sig_working;     %tank thickness, [ft]

rho_o = 1.4984;                 %liquid n2o density [sl/ft3]

V_o = m_o / rho_o;
V_caps_empty = 4/3 * pi * ((d_tank_o - 2*t_tank)/2)^3; %volume inside caps [ft3]
V_cyl_empty = V_o - V_caps_empty;       %vol inside cylinder [ft3]
L_cyl = V_cyl_empty / (pi / 4 * (d_tank_o - 2*t_tank)^2); %cyl length [ft]
L_tank = L_cyl + d_tank_o;            % tank length = cyl + 2 caps [ft]

% calc tank mass
V_caps_material = 4/3 * pi * ((d_tank_o/2)^3 - ((d_tank_o - 2*t_tank)/2)^3);
V_cyl_material = L_cyl * pi / 4 * (d_tank_o^2 - (d_tank_o - 2*t_tank)^2);
V_tank_material = V_caps_material + V_cyl_material;
m_tank = rho_Al * V_tank_material;      %tank mass [sl]
    

%outputs
t_tank = t_tank * 12;           %thickness tank [in]
m_tank = m_tank *32.2;          %tank mass [lbm]
end