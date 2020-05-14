function [rocket] = get_Buckling(rocket)
    % TODO: accout for combined load on ox tank
    % calculates neccessary wall thickness for body tube under axial load
    % assumed load is avg thrust plus drag at max V(but sea level air density. not accurate, but is an overestimate for safety)
    % Change j at start of script to switch to load bearing tanks (0,1) = (skin,skin+tank)

j=0; %to retain original function while modifying

%% Load Bearing Skin
if(j==0)
%note: freedom units (psf, lb/ft3)
d_body = rocket.geo.body.D; % [ft]
Yield = rocket.geo.body.material.yield; % [psf]
E = rocket.geo.body.material.young; % [psf]
ro = d_body/2; % [ft]
T = rocket.prop.T_avg;
rho = rocket.geo.body.material.density; % [lb/ft3]

n = 4; %number for end conditions (both fixed = 4)

% Stress Calc and Wall thickness solver
tWall = 0; %starting wall thickness [ft]
minFOS = 0;
reqFOS = 1.5;
i = 1; %count variable to avoid infinite
while (minFOS < reqFOS)
    tWall = tWall + 0.005/12; % [ft]
    aC = pi * (ro^2 - (ro-tWall)^2); % [ft^2]

    % Buckling length
    I = pi/4 * (ro^4 - (ro-tWall)^4);  % [ft^4]
    maxL = sqrt((n*pi^2*E*I) / (Yield*aC/minFOS)); % [ft]

    Tstress = T/aC; % [lb/ft2]
    Wstress = rho*maxL; % [lb/ft2]
    stress = Tstress + Wstress; % [lb/ft2]
    minFOS = Yield/stress;

    if(i == 10000)
        error('tWall not found')
    end

    i = i + 1;
end

tWall = max([tWall, 0.125/12]);
rocket.geo.body.l_buckle = maxL;
rocket.geo.body.t = tWall; % [ft]
end

%% Load Bearing Tanks
if(j==1)
    % assumes fixed ID(tank D) finds the OD and wall thickness
    % assumes that ox tank is loaded at the indicated FOS load
    % this would realistically translate to some surface area adapter for the ox tank to ensure that load is distributed as
    % intended within the rocket
    % TODO: Account for pressure caused hoop stress/combined load and relevant yield under such load to remove conservative FOS

    %Body tube data
    d_body = rocket.geo.ox_t.D; % [ft]
    l_body = rocket.geo.body.L;
    Yield = rocket.geo.body.material.yield; % [psf]
    E = rocket.geo.body.material.young; % [psf]
    ro = d_body/2; % [ft]
    rho = rocket.geo.body.material.density;
    n = 4; %number for end conditions (both fixed = 4)

    %Ox tank data
    d_tank = rocket.geo.ox_t.D;
    Yield_tank = rocket.geo.ox_t.material.yield;
    E_tank = rocket.geo.ox_t.material.young;
    ro_tank = d_tank/2; %[ft]
    ac_tank = pi*(ro_tank^2 - (ro_tank-rocket.geo.ox_t.t)^2);

    %Axial Load
    T = rocket.prop.T_avg;
    Swet_approx =  pi*ro^2;
    air = 0.0765; %density of air [lbm/ft^3] assuming sea level as a drag overestimate
    D = rocket.data.performance.D_max; %.5*air*rocket.data.performance.v_max^2*Swet_approx*rocket.aero.total.CD; %[lbm-ft/s^2]
    F = T+D;

    %% Stress Calc and Wall thickness solver
    minFOS = 1.5;
    %super conservative FOS for tank since atm I am not sure about accounting for combined loads(pressure hoop stress and axial compression considering they act in different directions, looks like it could be a material specific calculation)
    tankFOS = 2.5;

    % calculate force distributed to tank
    stress = Yield_tank/tankFOS + rho*l_body;
    F_tank = stress*ac_tank;

    % calculate body tube wall thickness and consequently, the OD of the body tube
    stress = Yield/minFOS + rho*l_body;
    aC = (F-F_tank)/stress; % [ft^2]
    ro = sqrt(aC/pi + ro_tank);
    tWall = ro-ro_tank;
    tWall = max([tWall, 0.125/12]);
    d_body = 2*ro;

    % Buckling length
    I = pi/4 * (ro^4 - (ro_tank)^4);  % [ft^4]
    I_tank = pi/4 * (ro_tank^4 - (ro_tank-rocket.geo.ox_t.t)^4);
    maxL = sqrt((n*pi^2*E*I) / (Yield*aC/minFOS)); % [ft]
    maxL_tank = sqrt((n*pi^2*E_tank*I_tank) / (Yield_tank*ac_tank/minFOS));


    rocket.geo.body.D = d_body;
    rocket.geo.body.l_buckle = maxL;
    rocket.geo.body.t = tWall; % [ft]
    rocket.geo.ox_t.l_buckle = maxL_tank;

end
end
