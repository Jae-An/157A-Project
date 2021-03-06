function rocket = rocket()
%Creates an instance of the plane class

rocket = struct('geo',[],'weight',[],'aero',[],'prop',[],'data',[]);

%% Rocket Geometry
rocket.geo = struct('nose',[],'body',[],'tail',[],'fins',[],'ox_t',[],'ox',[],'CC',[],'fuel',[],'payload',[],'misc',[]);
    rocket.geo.total = struct('L',[]);
    rocket.geo.nose = struct('D',[],'L',[],'t',[],'x',0,'A_p',[]);
    rocket.geo.body = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'K',1.5,'material',[]);
    rocket.geo.tail = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'R1',[],'R2',[]);
    rocket.geo.fins = struct('RC',[],'TC',[],'SS',[],'SL',[],'MC',[],'N',[],'t',[],'x',[]);
    
    rocket.geo.payload  = struct('D',0.25,'L',0.5,'x',[]);
    rocket.geo.ox_t     = struct('D',[],'L',[],'t',[],'x',[],'material',[]);
    rocket.geo.ox       = struct('D',[],'L',[],'x',[]);
    rocket.geo.CC       = struct('D',[],'L',[],'t',[],'x',[],'material',[]);
    rocket.geo.fuel     = struct('D_o',[],'D_i',[],'L',[]);
    rocket.geo.misc     = struct('L',[],'x',[]);
        % Contains 7 elements: main, recovery bulkhead, drogue, avionics, thrust
        % bulkhead, plumbing, nozzle (L and x will be arrays)
    
%% Weights and Inertias    
rocket.weight = struct('total',[],'nose',[],'body',[],'fins',[],'ox_t',[],'ox',[],'CC',[],'fuel',[],'misc',[]);
    rocket.weight.total = struct('W',[],'W_wet',[],'W_dry',[],'W_propellant',[],'CG',[]);
    rocket.weight.nose = struct('W',[],'CG',[]);
    rocket.weight.body = struct('W',[],'CG',[]);
    rocket.weight.tail = struct('W',[],'CG',[]);
    rocket.weight.fins = struct('W',[],'CG',[]);
    
    rocket.weight.payload = struct('W',[],'CG',[]);
    rocket.weight.press_t = struct('W',[],'CG',[]);
    rocket.weight.ox_t = struct('W',[],'CG',[]);
    rocket.weight.ox = struct('W',[],'CG',[]);  
    rocket.weight.CC = struct('W',[],'CG',[]);
    rocket.weight.fuel = struct('W',[],'CG',[]);  
    rocket.weight.misc = struct('W',[],'CG',[]);    

%% Aerodynamics
rocket.aero = struct('total',[],'nose',[],'body',[],'tail',[],'fins',[]);
    rocket.aero.total = struct('CD',[]); %
    rocket.aero.nose = struct('CD',[]);
    rocket.aero.body = struct('CD',[]);
    rocket.aero.tail = struct('CD',[]); 
    rocket.aero.fins = struct('CD',[]);  

%% Propulsion
rocket.prop = struct('T_avg',1000,'I',[],'Isp',[],'t_b',[],'OF',[],'A_e',[],'P_c',[],'h_opt',17000);

%% Performance, Aerodynamics, Stability
rocket.data = struct('requirements',[],'performance',[]);
    rocket.data.requirements = struct('v_ORS',100,'z_max',30000,'MOS',2); % Add stall angle?
    rocket.data.performance = struct('v_ORS',[],'z_max',[],'MOS',[],'Ma_max',[],'v_max',[],'a_max',[],'v_apogee',[],'D_max',[]);