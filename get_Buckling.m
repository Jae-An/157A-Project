function [rocket] = get_Buckling(rocket)
j=0; %to retain original function while modifying

if(j==0)
%% reference
%note: freedom units (psf, lb/ft3)
d_body = rocket.geo.body.D; % [ft]
Yield = rocket.geo.body.material.yield; % [psf]
E = rocket.geo.body.material.young; % [psf]
ro = d_body/2; % [ft]
T = rocket.prop.T_avg;
rho = rocket.geo.body.material.density; % [lb/ft3]

n = 4; %number for end conditions (both fixed = 4)

%% Stress Calc and Wall thickness solver
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

rocket.geo.body.l_buckle = maxL;
rocket.geo.body.t = tWall; % [ft]
end
%Load Bearing Tanks
if(j==1)
    %% reference
    %note: freedom units (psf, lb/ft3)
    %Body tube data
    d_body = rocket.geo.body.D; % [ft]
    Yield = rocket.geo.body.material.yield; % [psf]
    E = rocket.geo.body.material.young; % [psf]
    ro = d_body/2; % [ft]
    T = rocket.prop.T_avg;
    rho = rocket.geo.body.material.density; % [lb/ft3]
    n = 4; %number for end conditions (both fixed = 4)

    %Ox tank data
    d_tank = rocket.geo.ox_t.D;
    Yield_tank = rocket.geo.ox_t.material.yield;
    E_tank = rocket.geo.ox_t.material.young;
    ro_tank = d_tank/2; %[ft]




    %% Stress Calc and Wall thickness solver
    tWall = 0; %starting wall thickness [in]
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

        rocket.geo.body.l_buckle = maxL;
        rocket.geo.body.t = tWall; % [ft]
    end



end
end
