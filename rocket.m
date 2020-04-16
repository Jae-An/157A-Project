function rocket = rocket()
%Creates an instance of the plane class

rocket = struct('geo',[],'weight',[],'aero',[],'prop',[],'data',[]);

%% Rocket Geometry
rocket.geo = struct('nose',[],'body',[],'tail',[],'fins',[],'press_tank',[],'ox_tank',[],'CC',[],'misc',[]);
    rocket.geo.nose = struct('D',[],'L',[],'t',[],'x',0,'A_p',[],'A_w',[]);
    rocket.geo.body = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'K',1.5,'A_w',[]);
    rocket.geo.tail = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'R1',[],'R2',[],'A_w',[]);
    rocket.geo.fins = struct('RC',[],'TC',[],'SS',[],'SL',[],'MC',[],'N',[],'t',[],'x',[],'A_w',[]);
    rocket.geo.press_t = struct('D',[],'L',[],'t',[],'x',[]);
    rocket.geo.ox_t    = struct('D',[],'L',[],'t',[],'x',[]);
    rocket.geo.CC      = struct('D',[],'L',[],'t',[],'x',[]);
    rocket.geo.fuel    = struct('D_o',[],'D_i',[],'L',[]);

    rocket.geo.misc = struct('D',0.25,'L',0.5,'x',[]); % start with payload
    
%% Weights and Inertias    
rocket.weight = struct('total',[],'nose',[],'body',[],'tail',[],'fins',[],'press_tank',[],'ox_tank',[],'CC',[],'misc',[]);
    rocket.weight.total = struct('W',[],'W_wet',[],'W_dry',[],'W_propellant',[],'I',[],'CG',[]);
    rocket.weight.nose = struct('W',[],'I',[],'I_pt',[],'CG',[]);
    rocket.weight.body = struct('W',[],'I',[],'I_pt',[],'CG',[]);
    rocket.weight.tail = struct('W',[],'I',[],'I_pt',[],'CG',[]);
    rocket.weight.fins = struct('W',[],'I',[],'I_pt',[],'CG',[]);
    rocket.weight.press_t = struct('W',[],'I',[],'I_pt',[],'CG',[]);
    rocket.weight.ox_t = struct('W',[],'I',[],'I_pt',[],'CG',[]);
    rocket.weight.CC = struct('W',[],'I',[],'I_pt',[],'CG',[]);
    rocket.weight.fuel = struct('W',[],'W_i',[],'I',[],'I_pt',[],'CG',[],'CG_i',[]);
    rocket.weight.ox = struct('W',[],'W_i',[],'I',[],'I_pt',[],'CG',[],'CG_i',[]);    
    rocket.weight.misc = struct('W',[],'I',[],'I_pt',[],'CG',[]);    

%% Aerodynamics
rocket.aero = struct('total',[],'nose',[],'body',[],'tail',[],'fins',[]);
    rocket.aero.total = struct('CD',[]); %
    rocket.aero.nose = struct('CD',[]); %'aoa',[],'CNa',[],'L',[],'W',[],'m',[],'t',[],'CG',[],'A_p',[]);
    rocket.aero.body = struct('CD',[]); %'aoa',[],'CNa',[],'L',[],'W',[],'m',[],'t',[],'CG',[],'A_p',[],'K',1.5);
    rocket.aero.tail = struct('CD',[]); %'aoa',[],'CNa',[],'L',[],'W',[],'m',[],'t',[],'CG',[],'A_p',[],'R1',[],'R2',[]);
    rocket.aero.fins = struct('CD',[]); %'aoa',[],'CNa',[],'TC',[],'SS',[],'SL',[],'MC',[],'N',[],'W',[],'m',[],'t',[],'x',[],'CG',[]);
%     rocket.data.aero = struct('CL',[],'CL_alpha',[],'CD',[],'CD0',[],'CDi',[],'v_cruise',[],'LD',[],'CL_cruise',[]);
%     rocket.data.stability = struct('h_n',[],'is_stable',[],'yaw_is_stable',[],'stall',[],'alphas',[]);    

%% Propulsion
rocket.prop = struct('T_avg',[],'Isp',[],'t_b',[],'I',[],'OF',[],'A_e',[]);

%% Performance, Aerodynamics, Stability
rocket.data = struct('requirements',[],'performance',[]);
    rocket.data.requirements = struct('ORS',100,'z_max',30000); % Add stall angle?
    rocket.data.performance = struct('ORS',[],'z_max',[],'Ma_max',[],'a_max',[]);