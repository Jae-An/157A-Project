function [maxL,minFOS,tWall] = loadBody(rocket)
    %currently assuming that the material's own weight will not be significant
%TODO: get length where buckling becomes dominant failiure mode

%% reference
%note: freedom units
Yield = rocket.geo.body.material.yield;
E = rocket.geo.body.material.young;
ro = rocket.geo.body.D/2;
T = rocket.prop.T_avg;
rho = rocket.geo.body.material.density;

n = 4; %number for end conditions (both fixed = 4)

%% Stress Calc and Wall thickness solver
tWall = 0.035; %starting wall thickness [in]
aC = pi * (ro^2 - (ro-tWall)^2);
Tstress = rocket.prop.T_avg/aC; %longitudinal stress from thrust
Wstress = rho*g; % this assumes a unit long section of body tube
stress = Tstress + Wstress;
minFOS = Yield/stress;
i = 1; %count variable to avoid infinite
while (minFOS < 2 || minFOS > 2.3 && i < 10000)
    tWall = tWall + 0.005;
    aC = pi * (ro^2 - (ro-tWall)^2);
    Tstress = T/aC;
    Wstress = rho*g; % this assumes a unit long section of body tube
    stress = Tstress + Wstress;
    minFOS = Yield/stress;
end

if(i = 10000)
    maxL = 0;
    minFOS = 0;
    tWall = 0;
    return;
end
%% Buckling length
I = pi/64 * (rocket.geo.body.D^4 - (rocket.geo.body.D-tWall)^4);

maxL = sqrt((n*pi^2*E*I)/Yield);

end
