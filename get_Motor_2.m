% function [grain_d_i,grain_length, m_o, m_f, mdot,mdot_o, mdot_f_reqd, Thrust_jet, ...
%     Isp_new, nozzle_mass, A_e, A_t, epsilon, T_c, time_vec, Thrust_vec, Pe_vec] = ...
%     get_Motor_2(T, tb, Isp, P_c, OF ,d_o)
%OUTPUT grain inner diam [IN],grain length [FT], mass flow rates [SLUGS/S],
%Jet Thrust [LBF], Isp [S],nozzle mass [slugs], Exit Area [ft^2], throat
%area [ft^2],expansion ratio, exit pressure [PSF]
% ...propellant masses in SLUGS
% time_vec [S], thrust_vec [LBF]
% for an INPUT thrust [LBF], burn time {S], specific impulse [S]...
% chamber pressure [PSI], OF ratio, grain outer diam[IN]
% NOTE: returns zero inner grain diam if the input outer diam does not
% yield geometry with the required fuel mass flow rate 

% Inputs
T = 1000;
tb = 7.6;
Isp = 200;
P_c = 575;
OF = 6;
d_o = 6;

g = 32.2;

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
m_o = 0.86593;%mdot_o * tb;              %ox mass [sl]
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
%mdot_f = rho_f * A_b .* rdot;   % fuel mass flow rate [sl/s]
mdot_f = rho_f .* rdot * pi .*d_i .* L_grain;

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
P_e = 14.7*144; %assume fully expanded at 17000 ft (psf)
P_c = P_c * 144;
%V_rail = 100 ; %ft/s
%L_rail = 20; %ft

%design choices                   
mdot = mdot_tot; %[sl/s]                          %FROM Above

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
alpha = 18; %(degrees), Divergent cone half angle      %CHOSEN
alpha = alpha*pi/180; %Radians
cone_thickness = 1/100 *3.281; %ft %thickness of cone walls    %Guess
rho_nozzle = 1.8*100^3 /1000; %graphite density (kg/m^3)     %choice of material
rho_nozzle = rho_nozzle/515; %slug/ft^3
lamda = (1 + cos(alpha))/2;

%throat calculations
P_t = P_c*(1+ (k-1)/2)^(-k/(k-1)); %throat pressure
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

%% THRUST CURVE
% thrust curve for a static fire scenario, perf. expanded at sea level
%ASSUMPTIONS
% const. gamma = 1.25

%lookup table for pressure, temp, fluid quality as fn of rho, spec entropy
load('n2O_s_rho_props.mat');
%p_SR orig in bar
%t_SR orig in K
%q_SR dimensionless
%rho_vec kg/m3
%s_vec j/kg/K

p_SR = p_SR * 1e5 * 14.7 / 101325 * 144;      % pressure lookup array [psf]
rho_vec = rho_vec * 0.00194032;               % array of densities [sl/ft3]
s_vec = s_vec * 5.98;                  %spec. entropy array [ft lbf/sl R]


%%INPUTS

dt = 0.001;                   %time step [s]
t_tot = 10;                % sim time [s]
p_step = 10;              % pressure increment in pressure lookup vector

% tank and ox properties
rho_o_init = 1.4984;                 % ox density [sl/ft3]
%m_o = 1.0015;                   % ox mass for OF = 6, tb = 7.6 s
%mdot_o = 0.1318;                % ox flow rate [sl/s] hardcoded for now
V_o = m_o / rho_o_init;              % ox volume [ft3]
V_tank = V_o * 1.01;            % ox tank volume [+3% for ullage]
V_ullage = V_tank - V_o;        % ullage volume
T_tank = 536.4;                 %tank temp [deg R]

p_c = P_c;                      % chamber pressure [psf]
p_feed = 50 * 144;              % feed system delta p [psf]
p_inj = 0.2 * p_c;              % delta p injector [psf] (assume const)
%p_tank = p_c + p_inj + p_feed;             % tank pressure [psf]
p_tank = 745 * 144;

p_a = 14.7 * 144;               %ambient pressure [psf]



%fuel grain properties
%rho_htpb = 1.785;             %density cured htpb, [sl/ft3]
%rho_paraffin = 1.746;         %density paraffin, [sl/ft3]

%percent_htpb = .5;
%percent_paraffin = .5;

%rho_f = percent_htpb*rho_htpb + percent_paraffin*rho_paraffin; %[sl/ft3]

%mdot_tot = .1536;               % tot mass flow rate [sl/s]

mdot_f = mdot_f_reqd;                % fuel mdot [sl/s]
%m_f = 0.2167;             % fuel mass [sl] accounts for extra fuel
%V_f = m_f / rho_f;              %fuel volume [ft3]

d_o = d_o * 12;
d_i = grain_d_i;               %fuel grain diameters [in]
L_grain = grain_length;           % grain length [ft]

% regression rate
%a = 0.1146;                     %from paraffin/htpb/n2o research paper
%n = 0.5036;                     % for G_o in units of [lbm/in2/s]

% rocket properties
%p_e = 8.29   * 144;         % exit pressure, perf. expand at SL [psf]
%D_e = 2.46 /12;             % exit diam [ft]
%D_t = 0.88 / 12;            % throat diameter [ft]
%A_e = pi / 4 * D_e^2;       % exit area [ft2]
%A_t = pi / 4 * D_t^2;       % throat area [ft2]
%A_e = 0.0279;
%A_t = 0.0052;
%A_e = 0.047071732763238;    % nozzle for 575 psi pc
%A_t = 0.008694094039220;

eps = epsilon;
gam = k;                           %specific heat ratio (assume const.)

M_e = Ma_e;
pi_e = (1 + .5*(gam - 1)* M_e^2)^(-gam/(gam - 1)); %exit press ratio

d_elem = .05 / 12;           % injector hole diameter [ft]
n_inj = 50; %63                 % number of elements
A_elem = pi / 4 * d_elem^2; % area of each element [m2]
A_inj = n_inj * A_elem;     % injector area [m2]
Cd = 0.7;                   % discharge coeff
CdA = Cd * A_inj;

%%CALCULATIONS

%time loop
time_vec = 0:dt:t_tot;                  %total time vec, s
n_steps = length(time_vec);             %total number of time steps
step = 1;               %time step number

%preallocate arrays
Thrust_vec = zeros(n_steps,1);              %[lbf]
Mdot_vec = zeros(n_steps,1);                % tot mdot [sl/s]
OF_vec = zeros(n_steps,1);
Isp_vec = zeros(n_steps, 1);

Mdotf_vec = zeros(n_steps,1);               %fuel mdot [sl/s]
Mdoto_vec = zeros(n_steps,1);               %ox mdot [sl/s]
mo_vec = zeros(n_steps,1);                  %ox mass
mf_vec = zeros(n_steps,1);                  %fuel mass
Pt_vec = zeros(n_steps, 1);                 % tank press [psf]
Pc_vec = zeros(n_steps,1);                  % chamber press [psf]
Pe_vec = zeros(n_steps,1);                  % exit press [psf]

%assume perf. expand at 15000 ft
OFs = 1:1:13;                          %OF ratios
cstars = [3755, 4329, 4682, 4911, 5130, 5232, 5287, 5345, 5255, ... 
    5218,5141,5064,4992];
% char vel. assoc. w/ a certain of ratio 
% source: ProPep3



m_o_old = m_o;
m_f_old = m_f;
mdot_o_old = mdot_o;
mdot_f_old = mdot_f;
mdot_old = mdot_tot;
p_c_old = p_c;
p_t_old = p_tank;
p_e_old = P_e;
d_i_old = d_i;
A_cs_old = pi / 4 * (d_i_old )^2;    %cross-sec area[in2]
V_ullage_old = V_ullage;
T_tank_old = T_tank;
rho_o_old = rho_o_init;
V_o_old = V_o;

% initialize specific entropy & quality (based on pressure & init density)
s_old = 39 * 5.98;                         %[ft lbf / sl / R]
% spec. entropy obtained by goal seek method: find which entropy value
% returns a pressure of 745 psi when input in pressure lookup matrix
% p = interp2(s_vec, rho_vec, p_SR, s_old, rho_o_old, 'linear');
% initialize total entropy (ox mass * specific entropy)
s_tot = s_old * m_o_old;

% initialize quality 
quality = interp2(s_vec, rho_vec, q_SR, s_old, rho_o_old, 'linear');   

% Mdoto_vec(1) = mdot_o;
% Mdotf_vec(1) = mdot_f;
% Mdot_vec(1) = mdot_o + mdot_f;
% 
% Pc_vec(1) = p_c;
% Pt_vec(1) = p_tank;

while m_o_old >= 0  && p_t_old - p_c_old >= 0.20*p_c_old && quality <1 %run until ox is depleted
    t = time_vec(step);

    %calc coeff of jet thrust 
    I1 = 2*gam^2 / (gam - 1);
    I2 = 2 / (gam + 1);
    I3 = I2 ^ ((gam + 1) / (gam - 1));
    I4 = 1 - (pi_e) ^ ((gam - 1) / gam);
    c_tau = sqrt (I1 * I3 * I4);% + (p_e / p_c_old - p_a / p_c_old)*A_e / A_t;            %thrust coeff
    
    %estimate cstar
    OF = mdot_o_old / mdot_f_old;
    cstar = interp1(OFs,cstars, OF, 'linear', 5000);
    cstar = cstar * 0.8;           % 80% combustion efficiency
    %cstar = p_c_old * A_t / mdot_old;
    
    thrust = lamda*mdot_old * c_tau * cstar;    %thrust, lbf
    
    % calculations for next time step
    % new masses
    deltam_f = mdot_f_old * dt;
    m_f_new = m_f_old - deltam_f;
    deltam_o = mdot_o_old * dt;
    m_o_new = m_o_old - deltam_o;
    
    % find rate of change of entropy via specific entropy and ox mdot
    % get new total entropy
    s_tot = s_tot - deltam_o * s_old;
    % get new specific entropy
    s_new = s_tot / m_o_new;
    
    %calc new fuel mass flow rate
    G_o = mdot_o_old / A_cs_old * 32.2; %[lbm/s/in2]
    rdot = a * G_o.^n;              %regression rate [in/s]
    dr = rdot * dt;                 % change in grain diam [in]
    d_i_new = (d_i_old + 2*dr);     % new grain diam [in]
    
    A_cs_new = pi / 4 * d_i_new^2 / 144;    %new cross sect area [ft2]

    mdot_f_new = rho_f * rdot* pi *d_i_new * L_grain / 144;    %new fuel mass flow rate [sl/s]
    
    %calc new density based on new ox mass and vol
    V_o_new = m_o_new / rho_o_old;
    V_o_new = (V_o_new + V_o_old) / 2;      %take the average;
    rho_o_new = m_o_new / V_o_new;
    
    V_ullage_new = V_tank - V_o_new;
    
            % calc new tank pressure (USING LOOKUP TABLE AND NEW SPEC. ENTROPY
        % AND DENSITY)  (Scrap the below line)
    %p_t_new = p_t_old * (V_ullage_old / V_ullage_new)^1.27;
    p_t_new = interp2(s_vec, rho_vec, p_SR, s_new, rho_o_new, 'linear');
    
    %GET QUALITY
    quality = interp2(s_vec, rho_vec, q_SR, s_new, rho_o_new, 'linear');
    
    %calc tank temp assuming adiabatic cooling;
    T_tank_new = T_tank_old * (V_ullage_old / V_ullage_new)^(1.27 - 1);
    
    % calc new chamber pressure
    p_lookup = p_step:p_step:p_t_new;  %pressure lookup vec to solve for Pc
    mdot_out = p_lookup*A_t / cstar; %nozzle flow (flow out of chamber)
    dp = p_t_new - p_lookup;       %press differential
    mdot_in = CdA * sqrt(2 * rho_o_init * dp);   %injector flow (flow into chamber)
    delta_mdot = mdot_out - mdot_in;    % diff in nozzle and inj flowrates
    % look for chamber press in p_lookup that gives delta_mdot = mdot_f
    p_c_new = interp1(delta_mdot, p_lookup,mdot_f_new);
    
    %calc new exit pressure
    p_e_new = pi_e * p_c_new;
    
    %calc new ox mass flow rate and tot mass flow rate
    mdot_o_new = CdA * sqrt(2 * rho_o_init * (p_t_new - p_c_new));%% PUT rho_o_new LATER
    mdot_new = mdot_o_new + mdot_f_new;
    
    %log variables for plotting later
    Thrust_vec(step) = thrust;
    OF_vec(step) = OF;
    Isp_vec(step) = thrust / (mdot_o_old + mdot_f_old) / 32.2;
    Mdoto_vec(step) = mdot_o_old;
    Mdotf_vec(step) = mdot_f_old;
    Mdot_vec(step) = mdot_o_old + mdot_f_old;
    mo_vec(step) = m_o_old;
    mf_vec(step) = m_f_old;
    Pc_vec(step) = p_c_old;
    Pt_vec(step) = p_t_old;
    Pe_vec(step) = p_e_old;

    %update state variables
    m_f_old = m_f_new;
    m_o_old = m_o_new;
    mdot_f_old = mdot_f_new;
    mdot_o_old = mdot_o_new;
    mdot_old = mdot_f_old + mdot_o_old;
    p_t_old = p_t_new;
    p_c_old = p_c_new;
    p_e_old = p_e_new;
    d_i_old = d_i_new;
    A_cs_old = A_cs_new *144;       % update cross sect area to [in2]
    V_ullage_old = V_ullage_new;
    T_tank_old = T_tank_new;
    V_o_old = V_o_new;
    rho_o_old = rho_o_new;
    s_old = s_new;
    %update step number
    step = step + 1;
end

mo_vec(step:end) = m_o_old;
mf_vec(step:end) = m_f_old;

plot(time_vec, Thrust_vec);
xlabel('Time [s]');
ylabel('Thrust [lbf]');

save('Thrust_data','time_vec','Thrust_vec','Pe_vec','mo_vec','mf_vec')

%end