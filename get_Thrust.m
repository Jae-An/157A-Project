function [T_arr] = get_Thrust(time_vec, Thrust_vec, dt)
% Generates thrust array with uniform timesteps of dt using interpolation.
% T_data should be an xlsx file in the form 'FILENAME.xlsx'
% T_data should also start from t=0 and end at t=tb
% dt should be a power of 10

% Import data from file
t = time_vec(Thrust_vec>0);
T = Thrust_vec(Thrust_vec>0);

% Get tb
t_b = dt*floor(t(end)/dt); % no good way of rounding to nearest timestep (dt)

% Interpolate
tq = 0:dt:t_b;
T_arr = interp1(t,T,tq,'linear','extrap');

end