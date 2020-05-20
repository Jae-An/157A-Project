function [rocket] = get_Rocket_2(rocket)
% Fills rocket parameters, based on current design.

%% Geo
geo = rocket.geo;
    % Materials
    geo.body.material = mLibrary().Al.e; % Al 6061-T6
    geo.ox_t.material = mLibrary().Al.g; % Al 7050-T7651
    geo.CC.material   = mLibrary().Al.g; % Al 7050-T7651

    % Nose cone
    geo.nose.D = 6.25/12;                           % [ft] nose diameter
    FR = 5;                                         % [caliber] fineness ratio
    geo.nose.L = FR * geo.nose.D;                   % [ft] nose length
    geo.nose.t = 0.3/12;                            % [ft] nose thickness
    geo.nose.A_p = 0.5 * geo.nose.D * geo.nose.L;   % [ft2] nose projected area (~cone)
        
    % Body
    geo.body.D = geo.nose.D;                % [ft]
    AR = 20 - FR;                           % [caliber] body aspect ratio
    geo.body.L = AR * geo.body.D;           % [ft] body length
    geo.body.x = geo.nose.L;                % [ft] body starting position (wrt nose tip)
    geo.body.A_p = geo.body.D * geo.body.L; % [ft2] body projected area
        
    % Fins
    geo.fins.N = 4;                                     % [fins] fin count
    geo.fins.RC = 7.5/12;                               % [ft] fin root chord
    geo.fins.TC = 3.75/12;                              % [ft] fin tip chord
    geo.fins.SS = 4.125/12;                             % [ft] fin semispan
    geo.fins.SL = 7.25/12;                              % [ft] fin sweep length
    geo.fins.MC = sqrt(geo.fins.SS^2 + ...
        (geo.fins.SL+geo.fins.TC/2 - geo.fins.RC/2)^2); % [ft] fin mid chord line
    geo.fins.t = 0.06*(geo.fins.RC+geo.fins.TC)/2;      % [ft] fin average airfoil thickness
    geo.fins.x = geo.nose.L + geo.body.L - geo.fins.RC; % [ft] fins starting position (near end of rocket)
    
    % Aeroshell total
    geo.total.L = geo.nose.L + geo.body.L;  % [ft] total aeroshell length
    
    rocket.geo = geo;
    rocket = get_Buckling(rocket);          % returns body thickness [ft]
    
    %% Prop
    % Prop sizing
    prop = rocket.prop;

    prop.T_avg = 1000;                  % [lbf] average thrust
    prop.Isp = 200;                     % [s] specific impulse (initial estimate)
    prop.I = 7600;                      % [lbf*s] total impulse
    prop.t_b = prop.I / prop.T_avg;     % [s] burn time
    prop.OF = 6;                        % [] OF ratio
    prop.P_c = 575*144;                 % [psf] chamber pressure
    prop.h_opt = 17000;                 % [ft] altitude for optimization

    rocket.prop = prop;
    rocket.geo.ox_t.D = geo.body.D - 2*rocket.geo.body.t;   % Tank diameter is body inner diameter
    rocket = get_Motor_new(rocket);    % Sizes rest of propulsion system (tanks, CC, fuel, etc.)
    
    % Miscellaneous + Finishing prop geo
    % Order: nose bulkhead, drogue, payload bulkhead, avionics, main, bulkhead,
    % plumbing, nozzle

    geo = rocket.geo;
    
    geo.misc.L(1:7) = [0.5/12, 3.5/12, 0.5/12, 5.375/12, 8.5/12, 0.5/12, 0.5]; % [ft]
    
    geo.CC.x = geo.misc.x(8) - geo.CC.L;                % CC ends at nozzle
    geo.fuel.x = geo.CC.x + 0.05*geo.CC.L;              % Fuel in the middle of CC (CC.L = 1.1*fuel.L)
    geo.ox_t.x = geo.CC.x - geo.misc.L(7) - geo.ox_t.L; % Ox tank above CC, plumbing
    geo.ox.x = geo.ox_t.x + geo.ox_t.t;                 % Ox starts at ox tank
    
    geo.misc.x(1) = geo.nose.L * (3.252/12 / geo.nose.D);       % Fit nose bulkhead as far up nose as possible (D of 3.252")
    geo.misc.x(2) = geo.misc.x(1) + geo.misc.L(1);              % Drogue after nose bulkhead
    geo.misc.x(3) = geo.nose.L + geo.nose.D - geo.misc.L(4) - geo.payload.L - geo.misc.L(3); % Payload bulkhead in nose
    geo.misc.x(4) =  geo.nose.L + geo.nose.D - geo.misc.L(4);   % Avionics in nose shoulder
    geo.misc.x(5) = geo.misc.x(4) + geo.misc.L(4);              % Main chute after avionics
    geo.misc.x(6) = geo.misc.x(5) + geo.misc.L(5);              % Recovery bulkhead after main
    geo.misc.x(7) = geo.ox_t.x + geo.ox_t.L;                    % Plumbing after ox tank
    
    % Payload
    geo.payload.x = geo.misc.x(3) + geo.misc.L(3); % Payload after payload bulkhead
    
    rocket.geo = geo;

%% Weight [lb], CG [ft]
W = rocket.weight;       
	% Nose (Cone)
    W.nose.W = 115.492 * (0.5*pi*geo.nose.L*geo.nose.D*geo.nose.t); % Assumed fiberglass, conical nosecone, no tip
    W.nose.CG = 2*geo.nose.L/3; % Assume conical
    
    % Body (Cylinder shell)
    W.body.W = geo.body.material.density * pi*((geo.body.D/2)^2 - (geo.body.D/2 - geo.body.t)^2) * geo.body.L;
    W.body.CG = geo.body.x + geo.body.L/2;

    % Fins
    V_fins = 4 * 5.685 / 12^3;          % [ft3] from CAD
    W.fins.W = V_fins * 111.24;         % [lb] Uses carbon fiber density
    W.fins.CG = geo.fins.x + 0.4885;    % [ft] from CAD

    % Payload (Cylinder)
    W.payload.W = 2.2;
    W.payload.CG = geo.payload.x + geo.payload.L/2;

    % Misc (Cylinders (diameter doesn't matter))
    W.misc.W(1:7) = [1.7, 0.44, 0.42, 2.3, 0.92, 1.7, 3]; % [lb] initial estimates
    W.misc.CG = geo.misc.x + geo.misc.L/2;
    
    % Propulsion
    W.ox_t.CG = geo.ox_t.x + geo.ox_t.L/2;
    W.fuel.CG = geo.fuel.x + geo.fuel.L/2;
    W.CC.CG = geo.CC.x + geo.CC.L/2;
    W.ox.CG = geo.ox.x + geo.ox.L/2;
    
    % TOTALS
    W.total.W_wet = W.nose.W + W.body.W + W.fins.W ...
                  + W.ox_t.W + W.CC.W + W.fuel.W + W.ox.W + sum(W.misc.W);
    W.total.W_dry = W.total.W_wet - W.total.W_propellant; 
    
    W_arr_dry = [W.nose.W, W.body.W, W.fins.W, W.payload.W, W.ox_t.W, W.CC.W, sum(W.misc.W)];    
    CG_arr_dry = [W.nose.CG, W.body.CG, W.fins.CG, W.payload.CG, W.ox_t.CG, W.CC.CG, sum(W.misc.W.*W.misc.CG)/sum(W.misc.W)]; 
    W.total.CG_dry = sum(W_arr_dry.*CG_arr_dry) / W.total.W_dry;
    W.total.CG_wet = (sum(W_arr_dry.*CG_arr_dry) + W.fuel.W*W.fuel.CG + W.ox.W*W.ox.CG) / W.total.W_wet;
    
    rocket.weight = W;

%% CP and Margin of Stability
body = rocket.geo.body;
fins = rocket.geo.fins;
weight = rocket.weight;

    % Component CNalpha [rad-1]
    CNa_n = 2;
    CNa_b = 4*body.L*deg2rad(12.2)/(pi*body.D); % 12.2 deg is max ORAOA
    CNa_f = (1 + (body.D/2)/(fins.SS+(body.D/2))) * (4*fins.N*(fins.SS/body.D)^2 / (1+sqrt(1+(2*fins.MC/(fins.RC+fins.TC))^2)));
    
    % Component CP [ft]
    X_n = 0.666*rocket.geo.nose.L; % assume cone
    X_b = body.x + body.L/2;
    X_f = fins.x ...
        + fins.SL * (fins.RC + 2*fins.TC)/(3*(fins.RC+fins.TC)) ...
        + (fins.RC + fins.TC - fins.RC*fins.TC/(fins.RC+fins.TC))/6; 
    
    % Total CP [ft]
    X_CP_aoa = (X_n*CNa_n + X_b*CNa_b + X_f*CNa_f) / (CNa_n + CNa_b + CNa_f);
    X_CP_0   = (X_n*CNa_n + X_f*CNa_f) / (CNa_n + CNa_f);
    
    % Margin of Stability
    rocket.data.performance.MOS_aoa_wet = (X_CP_aoa - weight.total.CG_wet) / body.D;
    rocket.data.performance.MOS_aoa_dry = (X_CP_aoa - weight.total.CG_dry) / body.D;
    rocket.data.performance.MOS_0_wet = (X_CP_0 - weight.total.CG_wet) / body.D;
    rocket.data.performance.MOS_0_dry = (X_CP_0 - weight.total.CG_dry) / body.D;
end
