function [t_cc, d_cc, L_cc, m_cc] = get_cc(p_c,d_tank_o,grain_length)
%INPUT: chamber pressure [PSI],tank outer diam [IN], grain length [FEET]
%OUTPUT: cc thickness INCHES, cc diam INCHES, cc length FEET,
%  ... mass cc POUNDS

sig_yield = 276e6 * 14.7 / 101325 * 144;          %tensile strength, [psf]
FOS = 1.5;                  %factor of safety, CHOSEN VALUE
sig_working = sig_yield / FOS;  %working stress, [psf]
rho_Al = 5.2389;              %density, [sl/ft3]

p_c = p_c * 144;            %chamber pressure [psf]

d_cc = d_tank_o     / 12;                % cc diam [ft]             
L_cc = grain_length * 1.1;  %10% longer than fuel grain to allow mixing [ft]
t_cc = d_cc / 2 * p_c / sig_working;   %assume cc is same material as tank
%[ft]

V_cc_material = L_cc * pi / 4 * (d_cc^2 - (d_cc - 2*t_cc)^2);
m_cc = rho_Al * V_cc_material;          %cc mass [sl]

%OUTPUTS
t_cc = t_cc * 12;           %cc thickness [in]
d_cc = d_cc * 12;           %cc diameter [in]
m_cc = m_cc * 32.2;         %cc mass [lbm]
end