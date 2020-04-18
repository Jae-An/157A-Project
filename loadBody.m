function [maxL,minFOS,tWall] = loadBody(rocket)
    %currently assuming that the material's own weight will not be significant
%TODO: get length where buckling becomes dominant failiure mode

%% reference
%note: freedom units
Yield = rocket.geo.body.material.yield;
E = rocket.geo.body.material.young;

n = 4; %number for end conditions (both fixed = 4)

%% Stress Calc and Wall thickness solver
tWall = 0.035; %starting wall thickness [in]
aC = pi() * (rocket.geo.body.D^2 - (rocket.geo.body.D-tWall)^2);
stress = rocket.prop.T_avg/aC;

minFOS = Yield/stress;
while (minFOS < 2)
    tWall = tWall + 0.005;
    aC = pi() * (rocket.geo.body.D^2 - (rocket.geo.body.D-tWall)^2);
    stress = rocket.prop.T_avg/aC;

    minFOS = Yield/stress;
end

%% Buckling length
I = pi()/64 * (rocket.geo.body.D^4 - (rocket.geo.body.D-tWall)^4);

maxL = sqrt((n*pi^2*E*I)/Yield);

end
