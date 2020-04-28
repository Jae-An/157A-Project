% Runs a bunch of rockets and saves rockets with good performance
% Units: FPS
tic
clear variables; close all;
fprintf('Simulation Started \n')

%% Begin study
% New rocket instance
newRocket = rocket();
newRocket = get_Rocket(newRocket);
newRocket = get_Stability(newRocket);
newRocket = get_Performance(newRocket);

% Check for pass
if ~set_PassFail(newRocket)
    error('rocket no work')
end

% Plot
toc