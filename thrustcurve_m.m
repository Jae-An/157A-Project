clear all; clc; close all;

% thrust curve for a static fire scenario, perf. expanded at sea level
%ASSUMPTIONS
%perf expand. at SL, tank pressure decr. 30 psi/s (Self Pressurizing Tank
%   Dynamics, Zimmerman et al.)
% const. gamma = 1.25

%%INPUTS

dt = 0.001;                   %time step [s]
t_tot = 10;                % sim time [s]

%launch site parameters
g = 32.2;                   % grav accel, ft/s2


% tank and ox properties
rho_o = 1.4984;                 % ox density [sl/ft3]
m_o = 1.01235;                   % ox mass for OF = 6, tb = 7.6 s
mdot_o = 0.1318;                % ox flow rate [sl/s] hardcoded for now
V_o = m_o / rho_o;              % ox volume [ft3]
V_tank = V_o * 1.03;            % ox tank volume [+3% for ullage]
V_ullage = V_tank - V_o;        % ullage volume
T_tank = 536.4;                 %tank temp [deg R]

p_c = 575 * 144;                      % chamber pressure [psf]
p_feed = 50 * 144;              % feed system delta p [psf]
p_inj = 0.2 * p_c;              % delta p injector [psf] (assume const)
p_tank = p_c + p_inj + p_feed;             % tank pressure [psf]

p_a = 14.7 * 144;               %ambient pressure [psf]



%fuel grain properties
rho_htpb = 1.785;             %density cured htpb, [sl/ft3]
rho_paraffin = 1.746;         %density paraffin, [sl/ft3]

percent_htpb = .5;
percent_paraffin = .5;

rho_f = percent_htpb*rho_htpb + percent_paraffin*rho_paraffin; %[sl/ft3]

mdot_tot = .1536;               % tot mass flow rate [sl/s]

mdot_f = 0.022;                % fuel mdot [sl/s]
m_f = 0.2167;             % fuel mass [sl] accounts for extra fuel
V_f = m_f / rho_f;              %fuel volume [ft3]

d_o = 6;
d_i = 5.0518;               %fuel grain diameters [in]
L_grain = 2.1472;           % grain length [ft]

% regression rate
a = 0.1146;                     %from paraffin/htpb/n2o research paper
n = 0.5036;                     % for G_o in units of [lbm/in2/s]

% rocket properties
[~,p_e,~] = get_Atmosphere(0,17000);         % exit pressure, perf. expand at SL [psf]
D_e = 2.46 /12;             % exit diam [ft]
D_t = 0.88 / 12;            % throat diameter [ft]
A_e = pi / 4 * D_e^2;       % exit area [ft2]
A_t = pi / 4 * D_t^2;       % throat area [ft2]
A_e = 0.0279;
A_t = 0.0052;
eps = A_e / A_t;

%%CALCULATIONS

%time loop
time_vec = [0:dt:t_tot].';                  %total time vec, s
n_steps = length(time_vec);             %total number of time steps
    
%preallocate arrays
Thrust_vec = zeros(n_steps,1);              %[lbf]
Mdot_vec = zeros(n_steps,1);                % tot mdot [sl/s]
OF_vec = zeros(n_steps,1);

Mdotf_vec = zeros(n_steps,1);               %fuel mdot [sl/s]
Mdoto_vec = zeros(n_steps,1);               %ox mdot [sl/s]
Pt_vec = zeros(n_steps, 1);                 % tank press [psf]
Pc_vec = zeros(n_steps,1);                  % chamber press [psf]

%assume perf. expand at ~15000 ft
gam = 1.25;                           %specific heat ratio (assume const.)
OFs = 1:1:13;                          %OF ratios
cstars = [3755, 4329, 4682, 4911, 5130, 5232, 5287, 5345, 5255, ... 
    5218,5141,5064,4992];
% char vel. assoc. w/ a certain of ratio 
% source: ProPep3

step = 1;               %time step number


m_o_old = m_o;
m_f_old = m_f;
mdot_o_old = mdot_o;
mdot_f_old = mdot_f;
mdot_old = mdot_tot;
p_c_old = p_c;
p_t_old = p_tank;
d_i_old = d_i;
A_cs_old = pi / 4 * (d_i_old )^2;    %cross-sec area[in2]
V_ullage_old = V_ullage;
T_tank_old = T_tank;

Mdoto_vec(1) = mdot_o;
Mdotf_vec(1) = mdot_f;
Mdot_vec(1) = mdot_o + mdot_f;

Pc_vec(1) = p_c;
Pt_vec(1) = p_tank;

while m_o_old >= 0  && p_c_old <= p_t_old   %run until ox is depleted
    t = time_vec(step);
    
    %calc coeff of thrust
    I1 = 2*gam^2 / (gam - 1);
    I2 = 2 / (gam + 1);
    I3 = I2 ^ ((gam + 1) / (gam - 1));
    I4 = 1 - (p_e / p_c_old) ^ ((gam - 1) / gam);
    c_tau = sqrt (I1 * I3 * I4) + (p_e / p_c_old - p_a / p_c_old)*A_e / A_t;            %thrust coeff
    
    %estimate cstar
    OF = mdot_o_old / mdot_f_old;
    cstar(step) = interp1(OFs,cstars, OF, 'linear', 5000);
    cstar(step) = cstar(step) * 0.8;           % 90% combustion efficiency
    %cstar = p_c_old * A_t / mdot_old;
    
    thrust = mdot_old * c_tau * cstar;    %thrust, lbf
    
    % calculations for next time step
    % new masses
    deltam_f = mdot_f_old * dt;
    m_f_new = m_f_old - deltam_f;
    deltam_o = mdot_o_old * dt;
    m_o_new = m_o_old - deltam_o;
    
    
    %calc fuel mass flow rate
    G_o = mdot_o_old / A_cs_old * 32.2; %[lbm/s/in2]
    rdot = a * G_o.^n;              %regression rate [in/s]
    dr = rdot * dt;                 % change in grain diam [in]
    d_i_new = (d_i_old + 2*dr);     % new grain diam [in]
    
    A_cs_new = pi / 4 * d_i_new^2 / 144;    %new cross sect area [ft2]
    A_b = A_cs_new * L_grain;       %burn area [ft2]
    
    mdot_f_new = rho_f * A_b * rdot;    %new fuel mass flow rate [sl/s]
    
    % calc new ox mass flow rate
    V_o = m_o_new / rho_o;
    V_ullage_new = V_tank - V_o;        
    %p_t_new = p_t_old * (V_ullage_old / V_ullage_new)^1.27;
        %calc tank temp assuming adiabatic cooling;
    T_tank_new = T_tank_old * (V_ullage_old / V_ullage_new)^(1.27 - 1);
        %calc press temp 
    %p_t_new = p_t_old * (T_tank_old / T_tank_new)^(1.27 / (1-1.27));
    p_t_new = p_t_old - (30*144)*dt;        %Self press; research paper
    %p_c_new = mdot_old * cstar / A_t;
    p_c_new = p_t_new - p_feed - p_inj;  % assume feed and inj press const.
    %p_inj =  p_t_new - p_c_new - p_feed;
    
    mdot_o_new = 0.65 * (pi / 4 * (3/8)^2/ 144) * sqrt(2*rho_o * (p_t_new - p_c_new));
    
    %log variables for plotting later
    Thrust_vec(step) = thrust;
    OF_vec(step) = OF;
    Mdoto_vec(step+1) = mdot_o_new;
    Mdotf_vec(step+1) = mdot_f_new;
    Mdot_vec(step+1) = mdot_o_new + mdot_f_new;
    Pc_vec(step+1) = p_c_new;
    Pt_vec(step+1) = p_t_new;
    %T(step+1) = T_tank_old;

    %update state variables
    m_f_old = m_f_new;
    m_o_old = m_o_new;
    mdot_f_old = mdot_f_new;
    mdot_o_old = mdot_o_new;
    mdot_old = mdot_f_old + mdot_o_old;
    p_t_old = p_t_new;
    p_c_old = p_c_new;
    d_i_old = d_i_new;
    A_cs_old = A_cs_new *144;       % update cross sect area to [in2]
    V_ullage_old = V_ullage_new;
    T_tank_old = T_tank_new;
    %update step number
    step = step + 1;
end
%Thrust_vec(Thrust_vec > 1030) = 973;
%save('Thrust_data','time_vec','Thrust_vec')
