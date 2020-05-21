function [CD, Mach] = get_CD_Data(t_i, t_b, v_mag, atm_T, CD_data)

% Extract data
% Column 1 - 5: Mach, AOA, CD, CD power off, CD power on
data_Mach = CD_data(:,1);
data_CD_off = CD_data(:,4);
data_CD_on = CD_data(:,5);

% Find Mach number,
gamma = 1.401;
R = 1716.49; % air gas constant (ft lbf/slug R)
atm_T = atm_T + 491.67; % (R)
v_Sound = sqrt(gamma*R*atm_T); % (ft/s)
Mach = v_mag/v_Sound;

% Interpolate CD
if t_i <= t_b
    CD = interp1(data_Mach, data_CD_on, Mach, 'linear', 'extrap');
else
    CD = interp1(data_Mach, data_CD_off, Mach, 'linear', 'extrap');
end

end




