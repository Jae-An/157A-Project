function [Cd, Mach] = get_CD_Data(velocity, Temp, Cd_Excel, Mach_Excel)

% Find Speed of sound
gamma = 1.4;
R = 286;
Temp = (Temp*5/9-32)+273;  %Convert to K from F
v_Sound = sqrt(gamma*R*Temp)*3.28084;   %Calculate speed of sound in SI then convert to english

% Find Mach Number
Mach = velocity/v_Sound;

Cd = interp1(Mach_Excel, Cd_Excel, Mach, 'linear', 'extrap');

end




