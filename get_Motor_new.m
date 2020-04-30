function [grain_d_i,grain_length, m_o, m_f, mdot,mdot_o, mdot_f_reqd, Thrust_jet, ...
    Isp_new, nozzle_mass, A_e, A_t, epsilon, T_c,u_e] = ...
    get_Motor_new(T, tb, Isp,P_c, OF ,d_o)
%OUTPUT grain inner diam [IN],grain length [FT], mass flow rates [SLUGS/S],
%Jet Thrust [LBF], Isp [S],nozzle mass [slugs], Exit Area [ft^2], throat
%area [ft^2],expansion ratio
% ...propellant masses in SLUGS
% for an INPUT thrust [LBF], burn time {S], specific impulse [S]...
% chamber pressure [PSI], OF ratio, grain outer diam[IN]
% NOTE: returns zero inner grain diam if the input outer diam does not
% yield geometry with the required fuel mass flow rate 
g = 32.2;

%define parameters
%T = 1000;               %sea level thrust, [lbf]
%tb = 7.6;               %burn time
%Isp = 202.2;           %specific impulse, [s]  from ProPep3 software

%p_c = 575;              %chamber pressure, [psi]
%p_c = p_c * 144;        % chamber pressure, [psf]

%OF = 6;                 % ox fuel ratio

%cstar_theo = 5232.68;   % theoretical characteristic velocity, [ft/s]
%comb_eff = .95;         % combustion efficiency
%cstar = cstar_theo * .95;  %char. vel., [ft/s]

P_c = P_c*144; %psf

%fuel density
rho_htpb = 1.785;             %density cured htpb, [sl/ft3]
rho_paraffin = 1.746;         %density paraffin, [sl/ft3]

percent_htpb = .5;
percent_paraffin = .5;

rho_f = percent_htpb*rho_htpb + percent_paraffin*rho_paraffin; %[sl/ft3]
rho_o = 1.4984;                 %liquid n2o density [sl/ft3]

%calculate mass flow rate
mdot_tot = T / Isp / g; %total mass flow rate, sl/s
mdot_o = OF / (OF + 1) * mdot_tot;      %ox mass flowrate [sl/s]
mdot_f_reqd = 1 / (OF + 1) * mdot_tot;       %fuel mdot [sl/s]

%calculate mass
m_o = mdot_o * tb;              %ox mass [sl]
m_f = mdot_f_reqd * tb;              %fuel mass [sl]
theo_sliver_frac = 0.298;       %theoretical sliver fraction for circ. port
m_f = m_f * (1+theo_sliver_frac);   % account for leftover unburnt fuel

V_f = m_f / rho_f;              %fuel volume [ft3]

%initialize array of diameters
d_i = (0.1:0.1:d_o-0.1) / 12 ;          % inner diameter [ft]
d_o = d_o       /12;     %outer diameter [ft]
t_init = 0.5 * (d_o - d_i);     %initial thickness of grain, [ft]

% calc grain geometries
L_grain = V_f ./ (pi/4 * (d_o^2 - d_i.^2)); %grain length, [ft]
A_cs = pi ./ 4 * d_i.^2;         % grain port area   [ft2]
A_b = pi * d_i .* L_grain;       %grain burn area  [ft2]

LD = L_grain ./ d_o;            %length-diam ratio of grain
%stanford hybrid paper calls for L/D >= 4, for proper propellant mixing (?)

%% regression rate
a = 0.1146;                     %from paraffin/htpb/n2o research paper
n = 0.5036;                     % for G_o in units of [lbm/in2/s]
%a = 0.026;
%n = 0.8076;


G_o = mdot_o ./ A_cs;           %oxidizer mass velocity [sl/ft2/s]
G_o = G_o * 32.2 / 144;         %ox mass vel. [lbm/in2/s]

rdot = a * G_o.^n;              %regression rate [in/s]
rdot = rdot / 12;               %reg rate [ft/s]
mdot_f = rho_f * A_b .* rdot;   % fuel mass flow rate [sl/s]

% plot(d_i * 12, mdot_f*32.2);
% xlabel('grain inner diameter [in]');
% ylabel('fuel mass flow rate [lbm/s]')
% yline(mdot_f_reqd*32.2);    % intersects w/ graph to give 
% % the inner diameter that yields the required fuel mass flow rate

% figure;
% plot(d_i * 12, LD);
% xlabel('grain inner diameter [in]');
% ylabel('length-diameter-ratio');
% yline(4);   %graph values above this line have L/D ratio above 4

grain_d_i = interp1(mdot_f,d_i, mdot_f_reqd,'linear',0);
grain_length = V_f ./ (pi/4 * (d_o^2 - grain_d_i.^2));
grain_d_i = grain_d_i * 12;

%NOZZLE CALCULATIONS
%constants
R = 8.314; %universal gas constant (J/mol*K)
P_a = 14.7*144; %ambient pressure at sea level (psf)
P_e = 5.313*10^4; %assume fully expanded at 17000 ft (Pa)
P_e = P_e*0.021; %(psf)
%V_rail = 100 ; %ft/s
%L_rail = 20; %ft

%design choices                   
mdot = mdot_tot; %[lbm/s]                          %FROM Above

%combustion chamber parameters
R_specific = R/0.025427;   %[J/kg K]                     % from proprep3 program
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
lamda = (1 + cos(alpha))/2;

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

%Exit velocity
u_e = sqrt(2*k*R_specific*T_c/(k-1) *(1 - (P_e/P_c)^((k-1)/k))); %exit velocity (m/s)
%textbook value of u_e for HTPB/n20 at 1000 psi P_c = 1625 m/s

%thrust
Thrust_launch = lamda*(mdot*u_e + (P_e - P_a)*A_e); %thrust (lvf)
Thrust_jet = lamda*mdot*u_e;

%max launch weight
%Max_weight = Thrust_launch/(V_rail^2 /(2*L_rail*g) + 1); %in lbf
Isp_new = u_e/g; %specific impulse


end