function [rocket] = get_Performance(rocket)
    %% Initial values
    t_i = 0; t_f = 160; dt = 0.01;

    W_wet = rocket.weight.total.wet;
    W_dry = rocket.weight.total.dry;
    g = 32.174;

    T_avg = rocket.prop.T_avg;

    C_D = rocket.aero.total.CD;
    d_body = rocket.geo.body.D;
    A_c = pi*(d_body/2)^2;

    v_i = 0; z_i = 0;
    z_0 = 0; % Initial altitude above sea level

    %% Fill Initial Conditions
        t = t_i;

        W = W_wet;
        m = W / g;
        T = T_avg;

        v = v_i;
        x = z_i;

    %% Simulate

    while t <= t_f && v >= 0
        %% Current State Calculations
        % Atmosphere
        [atm_rho, ~, ~] = get_Atmosphere(x,z_0);

        % Drag force
        D = -sign(v) * (0.5*atm_rho*v^2) * C_D * A_c;

        % Newton's 2nd Law
        F = T - W + D;
        a = F / m;

        %% Calculate Future State
            if t+dt <= t_b
                W = W + W_dot*dt;
                m = W / g;
                T = T_avg;
            else
                W = W_dry;
                m = abs(W) / g;
                T = 0;
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
    