function [rocket] = get_Performance(rocket)
    %% Initial values
    t_i = 0; t_f = 160; dt = 0.01;

    W_wet = rocket.weight.total.W_wet;
    W_dry = rocket.weight.total.W_dry;
    g = 32.174;

    T_avg = rocket.prop.T_avg;
    Isp = rocket.prop.Isp;
    A_e = rocket.prop.A_e;
    
    d_body = rocket.geo.body.D;
    A_c = pi*(d_body/2)^2;

    v_i = 0; z_i = 0;
    z_0 = 0; % Initial altitude above sea level

    %% Fill Initial Conditions
        t = t_i;

        W = W_wet;
        m = W / g;

        v = v_i;
        x = z_i;

        [atm_rho_i, atm_P_i, atm_T_i] = get_Atmosphere(0,z_0);
        
        ORS = 0;
        
    %% Simulate

    while t <= t_f && v >= 0 && isreal(x)
        %% Current State Calculations
        % Atmosphere
        [atm_rho, atm_P, atm_T] = get_Atmosphere(x,z_0);

        % Thrust
        if W > W_dry
            T = T_avg + (atm_P_i-atm_P)*A_e;
        else
            T = 0;
        end
        
        % Weight
        W_dot = T / Isp;
        
        % Drag
        if (v ~= 0)
            C_D = get_CD(rocket, x+z_0, abs(v));
            q = 0.5 * atm_rho * v^2;
            D = -sign(v) * q * C_D * A_c;
        else
            D = 0;
        end
        
        % Newton's 2nd Law
        F = T - W + D;
        a = F / m;

        %% Calculate Future State
            if W - W_dot*dt > W_dry
                W = W + W_dot*dt;
                m = W / g;
            else
                W = W_dry;
                m = W / g;
            end

            % Euler method
            v = v + a*dt;
            x = x + v*dt;

            % Off-rail speed
            if x <= 20
                ORS = v;
            end
            
            % Move time forward
            t = t + dt;
    end
    
    %% Save data
    rocket.performance.z_max = x;
    rocket.performance.ORS = ORS;
    