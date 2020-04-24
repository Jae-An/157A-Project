function [CD, M] = get_CD(rocket, h, v)
% Note: v should be abs(v)
    %% General Parameters
        % Mach number
        if h <= 37000
            v_sound = -0.004*h + 1116.45;
        elseif 37000 < h && h <= 64000
            v_sound = 968.08;
        else
            v_sound = 0.0007*h + 924.99;
        end
        M = v/v_sound;
        
        % Kinematic viscosity
        if h <= 15000
            a = 0.00002503;
            b = 0;
        elseif 15000 < h && h <= 30000
            a = 0.00002760;
            b = -0.03417;
        else
            a = 0.00004664;
            b = -0.6882;
        end
        nu = 0.000157*exp(a*h + b);
        
        % Total length
        L = rocket.geo.total.L;
        
        % General surface roughness
        K = 0.00025;
        
    %% Friction Drag
        Cd_f = get_Cdf(rocket, v, M, nu, L, K);
        
    %% Base Drag
    L_b = rocket.geo.body.L;
    d = rocket.geo.body.D;
    K_b = 0.0274*atan(L_b/d + 0.0116);
    n = 3.6542 * (L_b/d)^-0.2733;
    
    if M <= 0.6
        Cd_b = K_b * (1)^n / sqrt(Cd_f); % fix this when we add boattail        
    else
        if 0.6 < M && M <= 1.0
            f_b = 1 + 215.8*(M - 0.6)^6;
        elseif 1.0 < M && M <= 2.0
            f_b = 2.0881*(M-1)^3 - 3.7938*(M-1)^2 + 1.4618*(M-1) + 1.883917;
        else
            f_b = 0.297*(M-2)^3 - 0.7937*(M-2)^2 + 0.1115*(M-2) + 1.64006;
        end
        Cd_b = f_b * (K_b * (1)^n / sqrt(get_Cdf(rocket, v, 0.6, nu, L, K)));
    end
        
    %% Wave Drag
    % Calculate Mach range
    L_n = rocket.geo.nose.L;
    M_D = -0.0156*(L_n/d)^2 + 0.136*(L_n/d) + 0.6817;
        
    if L_n/L < 0.2
        a = 2.4;
        b = -1.05;
    else
        a = -321.94*(L_n/L)^2 + 264.07*(L_n/L) - 36.348;
        b = 19.634 *(L_n/L)^2 - 18.369*(L_n/L) + 1.7434;
    end
    M_F = a*(L/d)^b + 1.0275;
    
    % Calculate Cd_w_max
    c = 50.676*(L_n/L)^2 - 51.734*(L_n/L) + 15.642;
    g = -2.2538*(L_n/L)^2 + 1.3108*(L_n/L) - 1.7344;
    if L/d >= 6
        Cd_w_max = c*(L/d)^g;
    else
        Cd_w_max = c*(6)^g;
    end
    
    % Calculate Cd_w
    if M < M_D
        Cd_w = 0;
    elseif M_D <= M && M < M_F
        x = (M - M_D) / (M_F - M_D);
        F = -8.3474*x^5 + 24.543*x^4 - 24.946*x^3 + 8.6321*x^2 + 1.1195*x;
        Cd_w = Cd_w_max * F;
    elseif M_F <= M
        Cd_w = Cd_w_max;
    end
    
    CD = Cd_f + Cd_b + Cd_w;
end