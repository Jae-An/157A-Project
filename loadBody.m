function [rocket] = get_Buckling(rocket)
    %currently assuming that the material's own weight will not be significant
%TODO: get length where buckling becomes dominant failure mode

%% reference
%note: freedom units (psi, lb/ft3)
d_body = rocket.geo.body.D * 12; % [in]
Yield = rocket.geo.body.material.yield; % [psi]
E = rocket.geo.body.material.young; % [psi]
ro = d_body/2;
T = rocket.prop.T_avg;
rho = rocket.geo.body.material.density / (12^3); % [lb/in3]

n = 4; %number for end conditions (both fixed = 4)

%% Stress Calc and Wall thickness solver
tWall = 0.035; %starting wall thickness [in]
aC = pi * (ro^2 - (ro-tWall)^2);
Tstress = rocket.prop.T_avg/aC; %longitudinal stress from thrust
Wstress = rho*g; % this assumes a unit long section of body tube
stress = Tstress + Wstress;
minFOS = Yield/stress;
reqFOS = 1.5;
i = 1; %count variable to avoid infinite
while (minFOS < reqFOS)% || minFOS > 2.3 && i < 10000)
    tWall = tWall + 0.005;
    aC = pi * (ro^2 - (ro-tWall)^2); % [in^2]
    
    % Buckling length
    I = pi/4 * (ro^4 - (ro-tWall)^4);
    maxL = sqrt((n*pi^2*E*I) / (Yield*aC/minFOS));
    
    Tstress = T/aC;
    Wstress = rho*maxL; 
    stress = Tstress + Wstress;
    minFOS = Yield/stress;
    
    i = i + 1;
end

if(i == 10000)
    maxL = 0;
    minFOS = 0;
    tWall = 0;
    return;
end

rocket.geo.body.

end
