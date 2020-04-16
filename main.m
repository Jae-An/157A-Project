% Problem 1 Mod Solution:
% One-Dimensional Flight Modeling
clear variables; close all;

%% Initial values
% Use values from Endurance
% Add initial/burn/final times, timestep, total timesteps
% Add wet/dry weight, gravity, average thrust
% Add initial conditions (a, v, z)

t_i = 0; t_f = 160; dt = 0.01; t_b = 16;

W_wet = -139.80;
W_dry = -91.53;
W_dot = (W_dry - W_wet)/t_b;
g = 32.174;

T_avg = 580;

C_D = 0.4;
d_body = 7/12;
A_c = pi*(d_body/2)^2;

v_i = 0; z_i = 0;
z_0 = 0; % Initial altitude above sea level

%% Initialize arrays
% We want to track:
% time
% weight, thrust, total force
% acceleration, velocity, position

m = 0;
W = 0;
T = 0;
D = 0;
F = 0;

atm_rho = 0;
atm_P = 0;
atm_T = 0;

a = 0;
v = 0;
x = 0;

%% Fill Initial Conditions
    t = t_i;

    W = W_wet;
    m = abs(W) / g;
    T = T_avg;
    
    v = v_i;
    x = z_i;
    
%% Simulate
% Follow pseudocode from slides

while t <= t_f && v >= 0
    %% Current State Calculations
    % Atmosphere
    [atm_rho, atm_P, atm_T] = get_Atmosphere(x,z_0);
    
    % Drag force
    D = -sign(v) * (0.5*atm_rho*v^2) * C_D * A_c;
    
    % Newton's 2nd Law
    F = T + W + D;
    a = F / m;
    
    %% Calculate Future State
        if t+dt <= t_b
            W = W + W_dot*dt;
            m = abs(W) / g;
            T = T_avg;
        else
            W = W_dry;
            m = abs(W) / g;
            T = 0;
        end

        % Euler method
        v = v + a*dt;
        x = x + v*dt;

        % Move time forward
        t = t+dt;
end