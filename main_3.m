% Runs rocket once to get performance.
% Units: FPS
tic
clear variables; close all;
fprintf('Simulation Started \n')

%% Begin study
% New rocket instance
newRocket = rocket();
newRocket = get_Rocket_2(newRocket);
newRocket = get_Performance(newRocket);

% Check for pass
if ~set_PassFail(newRocket)
    error('Rocket does not meet performance requirements.')
end


toc