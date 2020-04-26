function [rocket] = get_Performance(rocket)
% 2DOF sim. 
% Without rotational dynamics, this assumes that there is no wind
% but the rocket instantaneously changes to new flight path, pitch angle
% off the rail.

    %% Initial values
    length_rail = 20 - rocket.geo.total.L + rocket.weight.total.CG;
    theta_rail = deg2rad(-4);
    v_WG = [-22.005, 0].'; % 15 mph horizontal wind
    
    t_i = 0; t_f = 160; dt = 0.01;
    t_steps = (t_f-t_i)/dt + 1;
    t = linspace(t_i,t_f,t_steps);
    
    W_wet = rocket.weight.total.W_wet;
    W_dry = rocket.weight.total.W_dry;
    g = 32.174;

    T_avg = rocket.prop.T_avg;
    Isp = rocket.prop.Isp;
    A_e = rocket.prop.A_e;
    
    d_body = rocket.geo.body.D;
    A_c = pi*(d_body/2)^2;

    v_i = [0,0].';
    x_i = [0,0].';
    z_0 = 0; % Initial altitude above sea level
    
    p_i = theta_rail;
    
    %% Array initialization
    W = zeros(2,t_steps);
    T = zeros(2,t_steps);
        T_mag = zeros(1,t_steps);
    D = zeros(2,t_steps);
        C_D = zeros(1,t_steps);
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

    T_mag(1) = T_avg;
    
    v(:,1) = v_i;
    x(:,1) = x_i;

    p(1) = p_i;
    phi(1) = p_i;    
    off_rail = false;
    
    [atm_rho_i, atm_P_i, atm_T_i] = get_Atmosphere(0,z_0);

    %% Simulate
    i = 1;
    while t(i) <= t_f && v(2,i) >= 0 && isreal(x(:,i)) % apogee
        %% Current State Calculations
        % Atmosphere
        [atm_rho, atm_P, atm_T] = get_Atmosphere(x(i),z_0);
        
        % Thrust
        if W(2,i) > W_dry
            T_mag(i) = T_avg + (atm_P_i-atm_P)*A_e;
        else
            T_mag(i) = 0;
        end
        T(:,i) = T_mag(i) * [-sin(p(i)), cos(p(i))].';
        
        % Weight
        W_dot = T_mag(i) / Isp;
        
        % Drag
        q = 0.5 * atm_rho * rssq(v(:,i))^2;
        if rssq(v(:,i)) ~= 0
            [C_D(i), Mach(i)] = get_CD(rocket, x(i)+z_0, rssq(v(:,i)));
        else
            C_D(i) = 0;
            Mach(i) = 0;
        end
        D(:,i) = q * C_D(i) * A_c * [sin(phi(i)), -cos(phi(i))].';
            
        % Simple Rail Dynamics (No lift, thrust misalignment, drag is
        % parallel to rail
        if dist(i) <= length_rail
            A = [1, 1; rocket.geo.total.L-rocket.weight.total.CG, 0];
            B = [W(2,i)*sin(theta_rail); 0];
            R = sum(linsolve(A,B)) * [cos(theta_rail), sin(theta_rail)].';      
        else
            R = 0;
        end
        
        % Newton's 2nd Law
        F(:,i) = T(:,i) - W(:,i) + D(:,i) + R;
        a(:,i) = F(:,i) / m;

        %% Calculate Future State
        if W(2,i) - W_dot*dt > W_dry
            W(2,i+1) = W(2,i) - W_dot*dt;
            m = W(2,i+1) / g;
        else
            W(2,i+1) = W_dry;
            m = W(2,i+1) / g;
        end

        % Euler method
        v(:,i+1) = v(:,i) + a(:,i)*dt;
        x(:,i+1) = x(:,i) + v(:,i)*dt;
        dist(i+1) = dist(i) + abs(rssq(v(:,i))*dt);
        
        % Instantaneous rotation off rail
        if ~off_rail && dist(i+1)>length_rail
            v_RW = v(:,i+1) - v_WG;
            AOA = atan(-v_RW(1)/v_RW(2)) - theta_rail;
            RotMat = [cos(AOA), -sin(AOA); sin(AOA), cos(AOA)]; 
            v(:,i+1) = RotMat * v(:,i+1);
            p(i+1) = asin(-v(1,i+1)/rssq(v(:,i+1)));
            off_rail = true;
        else
            p(i+1) = p(i);
        end
        
        % Flight Path angle calculation
        if v(2,i+1)>=0
            phi(i+1) = asin(-v(1,i+1)/rssq(v(:,i+1)));
        else
            phi(i+1) = pi - asin(-v(1,i+1)/rssq(v(:,i+1)));
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
    
    
    
    rocket.data.performance.z_max = max(x(2,:));
    rocket.data.performance.Ma_max = max(Mach);
    rocket.data.performance.v_ORS = v_ORS;
    rocket.data.performance.v_max = max(rssq(v));
    rocket.data.performance.v_apogee = 
    rocket.data.performance.a_max = max(rssq(a));
    