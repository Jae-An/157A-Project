% Problem 3 Solution:
% 3DOF Flight Modeling
clc; clear variables; close all;
load('Thrust_data.mat')

%% Initial values
    % Environment/Launch pad
    g = 32.174;
        % Assume constant w/altitude
    v_WG = [-22.005; 0];
    rail_L = 20;
    rail_phi = deg2rad(-4);
    
    % Time
    t_i = 0; t_f = 300; dt = 0.01;
    t = [t_i:dt:t_f];
    t_steps = length(t);
    
    % Thrust (and Thrust_data interp)%%%%%%%%
    [T_jet, t_b] = get_Thrust(time_vec, Thrust_vec, dt);
    P_e = interp1(time_vec, Pe_vec, [t_i:dt:t_b]);
    m_o = interp1(time_vec, mo_vec, [t_i:dt:t_b]);
        m_o(t_steps) = m_o(end);
    m_f = interp1(time_vec, mf_vec, [t_i:dt:t_b]);
        m_f(t_steps) = m_f(end);
    nozzle_exit_d = (1.413*2) /12;
    nozzle_exit_A = 0.047071732763238;%pi*(nozzle_exit_d/2)^2;

    % Drag
    CD_data = readmatrix('CD Test.csv','Range','A2:E501');
    main_alt = 700;
    
    % Initial State
    v_i = [0; 0];
    x_i = [0; 0];
    z_0 = 0; % Initial altitude above sea level

    p_i = rail_phi; % Initial angle of rocket (and thus thrust), measured from vertical

%% Rocket geometry
% Use values from Endurance
    
% Dimensions [in]
    % Nose cone
    nose_d = 6.25 /12;
    nose_L = nose_d * 5;
    nose_AP = 0.5*nose_L*nose_d;
        % Modeled as a conical shape
    nose_CNa = 2;
    nose_CP = 0.666*nose_L;
    
    % Body
    body_d = nose_d;
    body_A = pi*(body_d/2)^2;
    body_L = body_d * 15;
    body_K = 1.5;
        % Factor used in Barrowman equations
    body_AP = body_L*body_d;
    body_X = nose_L;
        % Distance from nose tip to component
    body_CNa = 0;
        % No contribution for 0 AOA
    body_CP = nose_L + body_L/2;

    % Fins (Clipped Delta)
    fin_RC = 7 /12;
    fin_TC = 3.5 /12;
    fin_SS = 3.8 /12;
    fin_SL = 6.75 /12;
    fin_MC = sqrt(((fin_RC-fin_TC)/2)^2 + fin_SS^2);
        % Distance between root and tip mid-chords
    fin_N = 4;
        % Number of fins
    fin_X = nose_L + body_L - fin_RC;
        % Fin ends at end of rocket
    fin_CNa = (1 + (body_d/2)/(fin_SS+body_d/2)) ...
            * (4*fin_N*(fin_SS/nose_d)^2) ...
            / (1 + sqrt(1 + (2*fin_MC/(fin_RC+fin_TC))^2));
    fin_CP = fin_X + (fin_SL/3)*(fin_RC + 2*fin_TC)/(fin_RC + fin_TC) ...
           + (1/6)*(fin_RC + fin_TC - fin_RC*fin_TC/(fin_RC+fin_TC)); 
    
    % Fuel
    fuel_do = 6 /12;
    fuel_di_i = 5.0518 /12;
    fuel_L = 2.1472 /12;
    fuel_X = 90.158 /12;
    fuel_CG = fuel_X + fuel_L/2;
    fuel_rho = 1.7655;
    
    % Ox Tank
    oxt_d = 6 /12;
    oxt_t = 0.083 /12;
    oxt_X = 3.7375;
    oxt_L = (33.468/12) + oxt_d;
    % Oxidizer
    ox_d = oxt_d - 2*oxt_t;
    ox_A = pi*(ox_d/2)^2;
    ox_X_i = oxt_X + ox_d/6; % disregard initial ullage, treat as cylinder with same centroid, volume
    ox_L_i = (2*ox_d/3) + (oxt_L - oxt_d);
    ox_CG_i = ox_X_i + ox_L_i/2;
    ox_rho = 1.4984;
    
% Dry quantities
    W_dry = 81.5;  
    CG_dry = 74.99 /12;
    MOI_dry = 868.5;

% AP, CN, and CP
    i_AP = [nose_AP; body_AP];
    i_CNa = [nose_CNa; body_CNa; fin_CNa];
    i_CP = [nose_CP; body_CP; fin_CP];
        i_CP(2) = sum(i_CP(1:2).*i_AP) / sum(i_AP);
        
%% Initialize arrays
    % Complete mass properties (propellant mass known from thrust curve
    % script
    m = (W_dry/g) + m_f + m_o;      % [slug]    Mass
    W_mag = m * g;
    W = [zeros(1,t_steps);W_mag];  % [lb]      Weight vector
        W_f = m_f * g;
        W_o = m_o * g;
    
    % Oxidizer, CG
    ox_L = (m_o/ox_rho) / ox_A;
    ox_X = ox_X_i + ox_L_i - ox_L;
    ox_CG = ox_X + ox_L/2;
    CG = (CG_dry*W_dry + fuel_CG.*W_f + ox_CG.*W_o) ./ W_mag;
    
    % Fuel Inner Diameter
    fuel_A = (m_f/fuel_rho) / fuel_L;
    fuel_di = sqrt(fuel_do^2 - 4*fuel_A/pi);
    
    % MOI
    fuel_MOI_self = get_MOI_Cylinder(W_f,fuel_do,fuel_di,fuel_L);
    fuel_MOI = fuel_MOI_self + W_f.*(fuel_CG-CG).^2;
    ox_MOI_self = get_MOI_Cylinder(W_o,ox_d,0,ox_L);
    ox_MOI = ox_MOI_self + W_o.*(ox_CG-CG).^2;
    MOI = MOI_dry + fuel_MOI + ox_MOI;    
        
    T = zeros(2,t_steps);           % [lb]      Thrust vector
        T_mag = zeros(1,t_steps);   % [lb]      Thrust magnitude
        T_jet(t_steps) = 0;         % [lb]      Jet thrust (pads rest of jet thrust data array with zeros)
        P_e(t_steps) = 0;           % [lb/ft2]  Nozzle exit pressure
    N = zeros(2,t_steps);
    D = zeros(2,t_steps);           % [lb]      Drag vector
        CP = zeros(1,t_steps);      % [ft]      Center of pressure wrt nose tip
        CN = zeros(1,t_steps);      % []        Coefficient of normal force                
        CD = zeros(1,t_steps);      % []        Drag coefficient
    R = zeros(2,t_steps);           % [lb]      Reaction force from ground vector
    F = zeros(2,t_steps);           % [lb]      Net force vector
    M = zeros(1,t_steps);           % [lb*ft]   Net moment
    
    atm_rho = zeros(1,t_steps);     % [slg/ft3] Atmospheric density
    atm_P = zeros(1,t_steps);       % [lb/ft2]  Atmospheric pressure
    atm_T = zeros(1,t_steps);       % [F]       Atmospheric temperature

    q = zeros(1,t_steps);           % [lb/ft2]  Dynamic pressure
    Mach = zeros(1,t_steps);        % []        Mach number

    a = zeros(2,t_steps);           % [ft/s2]   Acceleration vectoy
    v = zeros(2,t_steps);           % [ft/s]    Velocity vector
    x = zeros(2,t_steps);           % [ft]      Position vector
    dist = zeros(1,t_steps);        % [ft]      Total distance traveled
    
    aa = zeros(1,t_steps);          % [rad/s2]  Angular acceleration
    w = zeros(1,t_steps);           % [rad/s]   Angular velocity
    p = zeros(1,t_steps);           % [rad]     Angular position
        phi = zeros(1,t_steps);     % [rad]     Flight path angle
        aoa = zeros(1,t_steps);     % [rad]     Angle of attack
        aoa_i = zeros(3,t_steps);   % [rad]     Component angle of attack
        
%% Fill Initial Conditions 
% (only for "calculate future step" variables)
    v(:,1) = v_i;
        v_RW = v_i;
    x(:,1) = x_i;
    
    p(:) = p_i;
    phi(1) = p_i;
    aoa(1) = 0;
    
%% Simulate
for i = 1:t_steps-1
    %% Current State Calculations
    % Atmosphere
    [atm_rho(i), atm_P(i), atm_T(i)] = get_Atmosphere(x(2,i),z_0);
    [CD(i), Mach(i)] = get_CD_Data(t(i), t_b, rssq(v(:,i)), atm_T(i), CD_data);
    
    % Thrust
    if t(i) <= t_b
        T_mag(i) = T_jet(i) + (P_e(i)-atm_P(i))*nozzle_exit_A; % correct for pressure thrust
    end
    T(:,i) = T_mag(i) * [-sin(p(i)), cos(p(i))].';
    
    % Normal Force
    q(i) = 0.5*atm_rho(i)*rssq(v_RW)^2;
    i_CNa(2) = 4*body_K*sum(i_AP)*aoa_i(2,i) / (pi*body_d^2);
    CN(i) = sum(i_CNa.*aoa_i(:,i));
    if aoa(i) ~= 0
        CP(i) = sum(i_CP.*abs(i_CNa)) / sum(abs(i_CNa));
    end
%     if CP(i) < CG(i) && dist(i) > rail_L
%         disp(t(i))
%         error('Statically unstable!')  
%     end
    N_mag = q(i) * CN(i) * body_A;
    N(:,i) = N_mag * [cos(p(i)); sin(p(i))];
    
    % Drag (No parachutes since only going to apogee for problem 3)
%     if v(2,i)<0 && t(i)>t_b
%         if x(2,i) > main_alt 
%             D_product = drogue_A * drogue_CD; % rocket body drag relatively negligible
%         elseif x(2,i) <= main_alt
%             D_product = main_A * main_CD; % drogue and rocket body drag relatively negligible
%         end
%     else
%         D_product = CD(i) * body_A; 
%     end
    D_product = CD(i) * body_A; 
    D(:,i) = q(i) * D_product * [sin(phi(i)); -cos(phi(i))];
    
    % Simple Rail Model (assume no rail friction)
    if dist(i) <= rail_L
        A = [1, 1; nose_L+body_L-CG(i), 0];
        B = [W(2,i)*sin(rail_phi); 0];
        R = sum(linsolve(A,B)) * [cos(rail_phi), sin(rail_phi)].'; % [lb] reaction force from rail   
    else
        R = 0;
    end
    
    % Newton's 2nd Laws
    F(:,i) = T(:,i) - W(:,i) + N(:,i) + D(:,i) + R;
    if (v(2,i)<0 && t(i)>t_b) || dist(i)<=rail_L
        M(i) = 0;
    else
        M(i) = (CP(i) - CG(i)) * N_mag;    
    end
        % Ground Reaction (contact) Force (from ground)
        if x(2,i)<=0 && F(2,i)<0
            G = [0; -F(2,i)];
            F(:,i) = F(:,i) + G;
        end
    a(:,i) = F(:,i) / m(i);
    aa(i) = M(i) / MOI(i);
    
    %% Calculate Future State
    % Euler method
    v(:,i+1) = v(:,i) + a(:,i)*dt;
    x(:,i+1) = x(:,i) + v(:,i)*dt;
    dist(i+1) = dist(i) + rssq(v(:,i))*dt;
    w(i+1) = w(i) + aa(i)*dt;
    p(i+1) = p(i) + w(i)*dt;
    
    % Additional angle/velocity calcs
    % Add rail model, better wind model later
    if dist(i+1) > rail_L
        v_RW = v(:,i+1) - v_WG;
        v_i_w = [(i_CP.'-CG(i+1))*w(i+1)*cos(p(i+1))...
                ;(i_CP.'-CG(i+1))*w(i+1)*sin(p(i+1))];
    else
        v_RW = v(:,i+1);
        v_i_w = zeros(2,3);
    end
    
    if rssq(v(:,i+1)) ~= 0
        if v(2,i+1)>=0
            phi(i+1) = asin(-v_RW(1)/rssq(v_RW));
        else
            phi(i+1) = pi - asin(-v_RW(1)/rssq(v_RW));
        end
    end

    aoa(i+1) = phi(i+1) - p(i+1);
    
    v_iW = v_RW + v_i_w;
    for j = 1:3
        if v_iW(2)>=0
            aoa_i(j,i+1) = asin(-v_iW(1,j)/rssq(v_iW(:,j))) - p(i+1);
        else
            aoa_i(j,i+1) = pi - (asin(-v_iW(1,j)/rssq(v_iW(:,j))) - p(i+1));
        end
    end    
    
    %% End Condition
    if t(i) > t_b && v(2,i) <= 0 % Break at apogee
        break
    end
end

%% Results
% Find apogee, range, max velocity, max acceleration
% as well as the times at which they happen

[x_max, i_x_max] = max(x(1,:));
[z_max, i_z_max] = max(x(2,:));
[v_max, i_v_max] = max(rssq(v));
[a_max, i_a_max] = max(rssq(a));

t_x_max = t(i_x_max);
t_z_max = t(i_z_max);
t_v_max = t(i_v_max);
t_a_max = t(i_a_max);

v_ORS = rssq(v);
v_ORS = v_ORS(dist > rail_L);
if ~isempty(v_ORS)
    v_ORS = v_ORS(1);
else
    v_ORS = 0;
end

%% Plotting
% Plot altitude vs. time
figure(1)
plot(t,x(2,:))
xlabel('Time (s)')
ylabel('Altitude (ft)')

% Plot velocity vs. time
figure(2)
plot(t,v(2,:))
xlabel('Time (s)')
ylabel('Vertical Velocity (ft/s)')

% Plot acceleration vs. time
figure(3)
plot(t,a(2,:))
xlabel('Time (s)')
ylabel('Vertical Acceleration (ft/s^{2})')

% Plot trajectory
figure(4)
plot(x(1,:),x(2,:))
title('Flight Trajectory')
xlabel('Distance (ft)')
ylabel('Altitude (ft)')
axis equal