function [rocket] = get_Performance(rocket)
% 2DOF sim. 
% Without rotational dynamics, this assumes that there is no angle of
% attack but the rocket instantaneously changes to new flight path, pitch angle
% off the rail.

    %% Initial values
    % General
    length_rail = 20 - rocket.geo.total.L + rocket.weight.total.CG_wet;
    theta_rail = deg2rad(-4);
    v_WG = [-22.005, 0].'; % 15 mph horizontal wind wrt ground inertial frame
    
    % Time
    t_i = 0; t_f = 1000; dt = 0.01;
    t_steps = (t_f-t_i)/dt + 1;
    t = linspace(t_i,t_f,t_steps);
    
    % Weight
    W_wet = rocket.weight.total.W_wet;
    W_dry = rocket.weight.total.W_dry;
    g = 32.174;
    
    % Thrust
    load('Thrust_data.mat')
    T_data = get_Thrust(time_vec,Thrust_vec,dt);
    T_data(t_steps) = 0;
    Isp = rocket.prop.Isp;
    A_e = rocket.prop.A_e;
    t_b = rocket.prop.t_b;
    
    % Drag
    Mach_arr = readmatrix('CD Test.csv','Range','A1:A500'); % RASAero II Data
    CD_off_arr = readmatrix('CD Test.csv','Range','F1:F500');
    CD_on_arr = readmatrix('CD Test.csv','Range','G1:G500');
    CD_drogue = 0.97;
    CD_main = 0.97;
    d_body = rocket.geo.body.D;
    d_drogue = 4;
    d_main = 10;
    A_c = pi*(d_body/2)^2;
    A_drogue = pi*(d_drogue/2)^2;
    A_main = pi*(d_main/2)^2;
    main_alt = 700; % main deployment altitude
    
    % Initial State
    v_i = [0,0].';
    x_i = [0,0].';
    z_0 = 0;
    h_opt = 0;%rocket.prop.h_opt;
    
    p_i = theta_rail;
    
    %% Array initialization
    W = zeros(2,t_steps);
    T = zeros(2,t_steps);
        T_mag = zeros(1,t_steps);
    D = zeros(2,t_steps);
        C_D = zeros(1,t_steps);
        C_D_prime = zeros(1,t_steps);
        q = zeros(1,t_steps);
    F = zeros(2,t_steps);
    
    Mach = zeros(1,t_steps);
    
    a = zeros(2,t_steps);
    v = zeros(2,t_steps);
    x = zeros(2,t_steps);
    dist = zeros(1,t_steps);
    
    p = zeros(1,t_steps); % angular position
    phi = zeros(1,t_steps); % flight path angle
    
    %% Fill Initial Conditions
    W(2,1) = W_wet;
    m = W(2,1) / g;

    T_mag(1) = T_data(1);
    
    v(:,1) = v_i;
    x(:,1) = x_i;

    p(1) = p_i;
    phi(1) = p_i;    
    off_rail = false;
    
    [~, P_e, ~] = get_Atmosphere(0,h_opt);

    %% Simulate
    i = 1;
    while t(i) <= t_f && x(2,i) >= 0 && isreal(x(:,i)) % apogee
        %% Current State Calculations
        % Atmosphere
        [atm_rho, atm_P, atm_T] = get_Atmosphere(x(2,i),z_0);
        
        % Thrust
        if W(2,i) > W_dry
            T_mag(i) = T_data(i) + (P_e-atm_P)*A_e;
        else
            T_mag(i) = 0;
        end
        T(:,i) = T_mag(i) * [-sin(p(i)), cos(p(i))].';
        
        % Weight
        W_dot = T_mag(i) / Isp;
        
        % Drag
        q(i) = 0.5 * atm_rho * rssq(v(:,i))^2;
        if rssq(v(:,i)) ~= 0
            if t(i) <= t_b
                [C_D(i), Mach(i)] = get_CD_Data(rssq(v(:,i)),atm_T,CD_on_arr,Mach_arr);
                [C_D_prime(i), ~] = get_CD(rocket, x(2,i)+z_0, rssq(v_RW));
            else
                [C_D(i), Mach(i)] = get_CD_Data(rssq(v(:,i)),atm_T,CD_off_arr,Mach_arr);
                [C_D_prime(i), ~] = get_CD(rocket, x(2,i)+z_0, rssq(v_RW));
            end
            if v(2,i)<0 && t(i)>t_b
                if x(2,i) > main_alt 
                    D_product = A_drogue * CD_drogue; % rocket body drag relatively negligible
                elseif x(2,i) <= main_alt
                    D_product = A_main * CD_main; % drogue and rocket body drag relatively negligible
                end
            else
                D_product = A_c * C_D(i);
            end
        else
            C_D(i) = 0;
            Mach(i) = 0;
            D_product = 0;
        end
        D(:,i) = q(i) * D_product * [sin(phi(i)), -cos(phi(i))].';
            
        % Simple Rail Dynamics (No lift, thrust misalignment, drag is
        % parallel to rail
        if dist(i) <= length_rail
            A = [1, 1; rocket.geo.total.L-rocket.weight.total.CG_wet, 0];
            B = [W(2,i)*sin(theta_rail); 0];
            R = sum(linsolve(A,B)) * [cos(theta_rail), sin(theta_rail)].'; % [lb] reaction force from rail   
        else
            R = 0;
        end
        
        % Newton's 2nd Law
        F(:,i) = T(:,i) - W(:,i) + D(:,i) + R;
        a(:,i) = F(:,i) / m;

        %% Calculate Future State
        if t(i+1) <= t_b
            W(2,i+1) = W(2,i) - W_dot*dt;
            m = W(2,i+1) / g;
        else
            W(2,i+1) = W_dry;
            m = W(2,i+1) / g;
        end

        % Euler method
        v(:,i+1) = v(:,i) + a(:,i)*dt;
        v_RW = v(:,i+1) - v_WG;
        x(:,i+1) = x(:,i) + v(:,i)*dt;
        dist(i+1) = dist(i) + abs(rssq(v(:,i))*dt);
        
        % Instantaneous rotation off rail
        if ~off_rail && dist(i+1)>length_rail % about to leave rail         
            AOA = atan(-v_RW(1)/v_RW(2)) - theta_rail;
            RotMat = [cos(AOA), -sin(AOA); sin(AOA), cos(AOA)]; 
            v(:,i+1) = RotMat * v(:,i+1);
            v_RW = v(:,i+1) - v_WG;
            p(i+1) = asin(-v(1,i+1)/rssq(v(:,i+1)));
            off_rail = true;
        else
            p(i+1) = p(i);
        end
        
        % Flight Path angle calculation
        if v(2,i+1)>=0
            phi(i+1) = asin(-v_RW(1)/rssq(v_RW));
        else
            phi(i+1) = pi - asin(-v_RW(1)/rssq(v_RW));
        end
        
        % Move time forward
        i = i + 1;
    end

    %% Save data
    v_ORS = rssq(v);
    v_ORS = v_ORS(dist > length_rail);
    if ~isempty(v_ORS)
        v_ORS = v_ORS(1);
    else
        v_ORS = 0;
    end
    
    [rocket.data.performance.z_max, i_a] = max(x(2,:));
    rocket.data.performance.t_apogee = t(i_a);
    rocket.data.performance.Ma_max = max(Mach);
    rocket.data.performance.v_ORS = v_ORS;
    rocket.data.performance.v_max = max(rssq(v));
    rocket.data.performance.v_apogee = rssq(v(:,i_a));
    rocket.data.performance.a_max = max(rssq(a));
    rocket.data.performance.D_max = max(rssq(D));
    
    %% Plotting
    RAS = readmatrix('Flight Test.csv');
    OPR = readmatrix('OpenRocket G10Sim.csv');
    % Trajectory
    figure(1)
    plot(x(1,1:i)/5280,x(2,1:i) , RAS(:,24)/5280,RAS(:,23),'--' , OPR(:,5)/5280,OPR(:,2),':k')
    xlabel('Distance [mi]')
    ylabel('Altitude [ft]')
    legend('Simulated', 'RASAero Data', 'OpenRocket Data')
    title('Flight Trajectory')
    grid on
    
    % Altitude vs. Time
    figure(2)
    plot(t(1:i),x(2,1:i), RAS(:,1),RAS(:,23),'--', OPR(:,1),OPR(:,2),':k')
    xlabel('Time [s]')
    ylabel('Altitude [ft]')
    title('Altitude Comparison')
    legend('Simulated', 'RASAero Data', 'OpenRocket Data')
    grid on
    
    % Mach (velocity) vs. Time
    figure(3)
    plot(t(1:i),Mach(1:i), RAS(:,1),RAS(:,4),'--', OPR(:,1),OPR(:,6),':k')
    xlabel('Time [s]')
    ylabel('Mach number')
    title('Mach Number Comparison')
    legend('Simulated', 'RASAero Data', 'OpenRocket Data')
    grid on
    
    % Acceleration vs. Time
    figure(4)
    plot(t(1:i),a(2,1:i), RAS(:,1),RAS(:,16),'--', OPR(:,1),OPR(:,4),':k')
    xlabel('Time [s]')
    ylabel('Acceleration [ft/s^{2}]')
    title('Vertical Acceleration Comparison')
    legend('Simulated','RASAero Data','OpenRocket Data')
    grid on
    
    % Pitch Angle vs. Time
    [~,i_a_data] = max(RAS(:,23)); 
    figure(5)
    plot(t(1:i_a),rad2deg(p(1:i_a))+90 ,  RAS((1:i_a_data),1),RAS((1:i_a_data),21),'--', ...
         t(1:i_a),rad2deg(phi(1:i_a))+90, RAS(1:i_a_data,1),RAS(1:i_a_data,22),'--')
    legend('Pitch (Sim)','Pitch (RAS)', 'Flight Path Angle (Sim)','Flight Path Angle (RAS)')
    grid on
    
    
    