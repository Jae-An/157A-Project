function [rocket] = get_Rocket(rocket)
% Fills new rocket parameters, based on initial OpenRocket.
% Final version will have rand ranges around initial values


%% Geo
geo = rocket.geo;
    % Materials
    geo.body.material = mLibrary().Al.e;
    geo.ox_t.material = mLibrary().Al.g;
    geo.CC.material = mLibrary().Al.g;

    % rocket.geo.nose = struct('D',[],'L',[],'t',[],'x',0,'A_p',[]);
        geo.nose.NC = 1; % VKO
        geo.nose.D = 6/12; %(4 + rand*(6 - 4)) / 12; % [ft]
        FR = 5; % + rand*(5 - 3);
        geo.nose.L = FR * geo.nose.D;
        geo.nose.t = 0.3/12; % keep constant
        geo.nose.A_p = 0.5 * geo.nose.D * geo.nose.L;
    % rocket.geo.body = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'K',1.5);
        geo.body.D = geo.nose.D; % [ft]
        AR = 20 - FR; %((15 + rand*(20 - 15)) - FR);
        geo.body.L = AR * geo.body.D;
        geo.body.x = geo.nose.L;
        geo.body.A_p = geo.body.D * geo.body.L;
    % rocket.geo.tail = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'R1',[],'R2',[]);
%         geo.tail.D = ; no tail for now
%         geo.tail.L = ;
%         geo.tail.t = ;
%         geo.tail.x = ;
%         geo.tail.A_p = ;
%         geo.tail.R1 = ;
%         geo.tail.R2 = ;
    % rocket.geo.fins = struct('RC',[],'TC',[],'SS',[],'SL',[],'MC',[],'N',[],'t',[],'x',[]);
        geo.fins.N = 4; %round(rand + 3); % sets N to 3 or 4
        geo.fins.x = geo.nose.L + geo.body.L; % initial approximation
    % Aeroshell total
    geo.total.L = geo.nose.L + geo.body.L; 
    
    % Payload
    geo.payload.x = geo.nose.L * (geo.payload.D/geo.nose.D);
    
    rocket.geo = geo;
    rocket = get_Buckling(rocket); % geo.body.t
    
    %% Prop
    % Prop sizing
        prop = rocket.prop;

        prop.T_avg = 1000; % [lbf] 
        prop.Isp = 200;
        prop.I = 6500;%4000 + rand*(8000 - 4000); % will fix
        prop.t_b = prop.I / prop.T_avg;
        prop.OF = 6; % Assume not varied
        prop.P_c = 500*144; % [psf] Assumed for now
        prop.h_opt = 0;
        
        rocket.prop = prop;
        rocket.geo.ox_t.D = geo.body.D - 2*rocket.geo.body.t; % For now approximate as wide tanks
        rocket = get_Motor(rocket);
    
    % Miscellaneous + Finishing prop geo
    % (main, rec. bulkhead, drogue, avionics, thr bulkhead, plumbing)    
    geo = rocket.geo;
    geo.misc.L(1:6) = [0.75, 0.0625, 0.5, 0.5, 0.05, 0.25]; % ALL ESTIMATES
    geo.CC.x = geo.total.L - geo.CC.L;
    geo.fuel.x = geo.CC.x + 0.05*geo.CC.L; % fuel in middle of CC
    geo.ox_t.x = geo.CC.x - geo.misc.L(5) - geo.ox_t.L;
    geo.ox.x = geo.ox_t.x + 2*geo.ox_t.t;
    
    geo.misc.x = geo.misc.L;
    dx = 1/12; % [ft] Spacing between consecutive components
    geo.misc.x(1) = geo.payload.x + geo.payload.L + dx;
    for i = 2:4
        geo.misc.x(i) = geo.misc.x(i-1) + geo.misc.L(i-1) + dx;
    end
    geo.misc.x(5) = geo.ox_t.x - geo.misc.L(4);
    geo.misc.x(6) = geo.ox_t.x + geo.ox_t.L + geo.misc.L(5);
    

    
    rocket.geo = geo;



%% Weight
W = rocket.weight;       
	% Nose (Cone)
    W.nose.W = 115.492 * (0.5*pi*geo.nose.L*geo.nose.D*geo.nose.t); % Assumed fiberglass, conical nosecone, no tip
    W.nose.CG = 2*geo.nose.L/3;
    
    % Body (Cylinder shell)
    W.body.W = geo.body.material.density * 0.25*pi*geo.body.L*(2*geo.body.D*geo.body.t - geo.body.t^2);
    W.body.CG = geo.body.x + geo.body.L/2;
    
    % Tail (empty for now)
    
    % Fins (Pt mass)
    W.fins.W = 4; % Guess from Endurance, conservative estimate
    W.fins.CG = geo.fins.x; % At end of rocket, conservative approximation

    % Payload (Cylinder)
    W.payload.W = 2.2;
    W.payload.CG = geo.payload.x + geo.payload.L/2;

    % Misc (Cylinders (diameter doesn't matter))
    W.misc.W(1:6) = [2.3, 3, 1.5, 2.5, 1.5, 4];
    W.misc.CG = geo.misc.x + geo.misc.L/2;
    % main, rec. bulkhead, drogue, avionics, thr bulkhead, plumbing
    
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

end