function rocket = rocket()
%Creates an instance of the plane class

rocket = struct('geo',[],'weight',[],'aero',[],'prop',[],'data',[]);

%% Rocket Geometry
rocket.geo = struct('nose',[],'body',[],'tail',[],'fins',[],'press_tank',[],'ox_tank',[],'plumbing',[],'CC',[],'payload',[],'misc',[]);
    rocket.geo.total = struct('L',[]);
    rocket.geo.nose = struct('D',[],'L',[],'t',[],'x',0,'A_p',[]);
    rocket.geo.body = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'K',1.5);
    rocket.geo.tail = struct('D',[],'L',[],'t',[],'x',[],'A_p',[],'R1',[],'R2',[]);
    rocket.geo.fins = struct('RC',[],'TC',[],'SS',[],'SL',[],'MC',[],'N',[],'t',[],'x',[]);
    
    rocket.geo.payload  = struct('D',0.25,'L',0.5,'x',[]);
    rocket.geo.ox_t     = struct('D',[],'L',[],'t',[],'x',[]);
    rocket.geo.ox       = struct('D',[],'L',[],'x',[]);
    rocket.geo.CC       = struct('D',[],'L',[],'t',[],'x',[]);
    rocket.geo.fuel     = struct('D_o',[],'D_i',[],'L',[]);
    rocket.geo.misc     = struct('L',[],'x',[]);
        % Contains main, recovery bulkhead, drogue, avionics, thrust
        % bulkhead, plumbing (L and x will be arrays)
    
%% Weights and Inertias    
rocket.weight = struct('total',[],'nose',[],'body',[],'tail',[],'fins',[],'press_tank',[],'ox_tank',[],'CC',[],'misc',[]);
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
    rocket.aero.nose = struct('CD',[]); %'aoa',[],'CNa',[],'L',[],'W',[],'m',[],'t',[],'CG',[],'A_p',[]);
    rocket.aero.body = struct('CD',[]); %'aoa',[],'CNa',[],'L',[],'W',[],'m',[],'t',[],'CG',[],'A_p',[],'K',1.5);
    rocket.aero.tail = struct('CD',[]); %'aoa',[],'CNa',[],'L',[],'W',[],'m',[],'t',[],'CG',[],'A_p',[],'R1',[],'R2',[]);
    rocket.aero.fins = struct('CD',[]); %'aoa',[],'CNa',[],'TC',[],'SS',[],'SL',[],'MC',[],'N',[],'W',[],'m',[],'t',[],'x',[],'CG',[]);  

%% Propulsion
rocket.prop = struct('T_avg',[],'Isp',[],'t_b',[],'OF',[],'A_e',[]);

%% Performance, Aerodynamics, Stability
rocket.data = struct('requirements',[],'performance',[]);
    rocket.data.requirements = struct('ORS',100,'z_max',30000,'MOS',1.25); % Add stall angle?
    rocket.data.performance = struct('ORS',[],'z_max',[],'MOS',[],'Ma_max',[],'a_max',[]);