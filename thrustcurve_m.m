clear all; clc; close all;

% thrust curve for a static fire scenario, perf. expanded at sea level
%ASSUMPTIONS
%perf expand. at SL, tank pressure decr. 30 psi/s (Self Pressurizing Tank
%   Dynamics, Zimmerman et al.)
% const. gamma = 1.25

%%INPUTS

load('n2O_s_rho_props.mat');
%p_SR orig in bar
%t_SR orig in K
%q_SR dimensionless
%rho_vec kg/m3
%s_vec j/kg/K

p_SR = p_SR * 1e5 * 14.7 / 101325 * 144;      % pressure lookup array [psf]
rho_vec = rho_vec * 0.00194032;               % array of densities [sl/ft3]
s_vec = s_vec * 5.98;                  %spec. entropy array [ft lbf/sl R]


dt = 0.001;                   %time step [s]
t_tot = 10;                % sim time [s]
p_step = 10;              % pressure increment in pressure lookup vector


%launch site parameters
g = 32.2;                   % grav accel, ft/s2


% tank and ox properties
rho_o_init = 1.4984;                 % ox density [sl/ft3]
m_o = 1.011535048802129;                   % ox mass for OF = 6, tb = 7.6 s
mdot_o = 0.1318;                % ox flow rate [sl/s] hardcoded for now
V_o = m_o / rho_o_init;              % ox volume [ft3]
V_tank = V_o * 1.01;            % ox tank volume [+1% for ullage]
V_ullage = V_tank - V_o;        % ullage volume
T_tank = 536.4;                 %tank temp [deg R]

p_c = 575 * 144;                      % chamber pressure [psf]
p_feed = 50 * 144;              % feed system delta p [psf]
p_inj = 0.2 * p_c;              % delta p injector [psf] (assume const)
%p_tank = p_c + p_inj + p_feed;             % tank pressure [psf]
p_tank = 745 * 144;

p_a = 14.7 * 144;               %ambient pressure [psf]



%fuel grain properties
rho_htpb = 1.785;             %density cured htpb, [sl/ft3]
rho_paraffin = 1.746;         %density paraffin, [sl/ft3]

percent_htpb = .5;
percent_paraffin = .5;

rho_f = percent_htpb*rho_htpb + percent_paraffin*rho_paraffin; %[sl/ft3]

mdot_tot = .1536;               % tot mass flow rate [sl/s]

mdot_f = 0.022;                % fuel mdot [sl/s]
m_f = 0.218828748890861;             % fuel mass [sl] accounts for extra fuel
V_f = m_f / rho_f;              %fuel volume [ft3]

d_o = 6;
d_i = 5.0518;               %fuel grain diameters [in]
L_grain = 2.1472;           % grain length [ft]

% regression rate
a = 0.1146;                     %from paraffin/htpb/n2o research paper
n = 0.5036;                     % for G_o in units of [lbm/in2/s]

% rocket properties
%[~,p_e,~] = get_Atmosphere(0,17000);         % exit pressure, perf. expand at SL [psf]
D_e = 2.46 /12;             % exit diam [ft]
D_t = 0.88 / 12;            % throat diameter [ft]
A_e = pi / 4 * D_e^2;       % exit area [ft2]
A_t = pi / 4 * D_t^2;       % throat area [ft2]
A_e = 0.075557076441258;
A_t = 0.008781034979612;       % nozzle for 575 psi pc
eps = A_e / A_t;
gam = 1.25111;                 %specific heat ratio (assume const.)

M_e = nozzleMach(eps,gam,"sup");
pi_e = (1 + .5*(gam - 1)* M_e^2)^(-gam/(gam - 1)); %exit press ratio

%injector properties (SUBJECT TO CHANGE)
%injector data
d_elem = .05 / 12;           % injector hole diameter [ft]
n_inj = 50; %63                 % number of elements
A_elem = pi / 4 * d_elem^2; % area of each element [m2]
A_inj = n_inj * A_elem;     % injector area [m2]
Cd = 0.7;                   % discharge coeff
CdA = Cd * A_inj;


%%CALCULATIONS

%time loop
time_vec = [0:dt:t_tot].';                  %total time vec, s
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

%assume perf. expand at ~15000 ft
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
p_e_old = pi_e * p_c;
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

while m_o_old >= 0  && p_t_old - p_c_old >= 0.15*p_c_old && quality <1 %run until ox is depleted
    t = time_vec(step);
    
    %calc coeff of thrust
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
    
    thrust = mdot_old * c_tau * cstar;    %thrust, lbf
    
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

save('Thrust_data','time_vec','Thrust_vec','Pe_vec','mo_vec','mf_vec')
