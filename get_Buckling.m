function [rocket] = get_Buckling(rocket,j)
    %currently assuming that the material's own weight will not be significant
%TODO: get length where buckling becomes dominant failure mode

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
while (minFOS < reqFOS)% || minFOS > 2.3 && i < 10000)
    tWall = tWall + 0.005/12; % [ft]
    aC = pi * (ro^2 - (ro-tWall)^2); % [ft^2]

    % Buckling length
    I = pi/4 * (ro^4 - (ro-tWall)^4);  % [ft^4]
    maxL = sqrt((n*pi^2*E*I) / (Yield*aC/minFOS)); % [ft]

    Tstress = T/aC; % [lb/ft2]
    Wstress = rho*maxL; % [lb/ft2]
    stress = Tstress + Wstress; % [lb/ft2]
    minFOS = Yield/stress;

    i = i + 1;
end

if(i == 10000)
    error('tWall not found')
end

rocket.geo.body.t = tWall; % [ft]
end
if(j==1)
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
    while (minFOS < reqFOS)% || minFOS > 2.3 && i < 10000)
        tWall = tWall + 0.005/12; % [ft]
        aC = pi * (ro^2 - (ro-tWall)^2); % [ft^2]

        % Buckling length
        I = pi/4 * (ro^4 - (ro-tWall)^4);  % [ft^4]
        maxL = sqrt((n*pi^2*E*I) / (Yield*aC/minFOS)); % [ft]

        Tstress = T/aC; % [lb/ft2]
        Wstress = rho*maxL; % [lb/ft2]
        stress = Tstress + Wstress; % [lb/ft2]
        minFOS = Yield/stress;

        i = i + 1;
    end

    if(i == 10000)
        error('tWall not found')
    end

    rocket.geo.body.t = tWall; % [ft]

end
end
