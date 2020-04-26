% Runs a bunch of rockets and saves rockets with good performance
% Units: FPS
tic
clear variables; close all;
fprintf('Optimization Started \n')

%% Study parameters
g = 0; % passed rockets
b = 0; % failed rockets
n = 0; % total# rockets
G = 1000; % desired# passed rockets
results = struct();

%% Begin study
fprintf('Finding good rockets... \n')
while g < G
    % New rocket instance3
    newRocket = rocket();
    newRocket = get_Rocket(newRocket);
    newRocket = get_Stability(newRocket);
    newRocket = get_Performance(newRocket);
    
    % Check for pass
    if set_PassFail(newRocket)
        g = g + 1;
        results(g).rocket = newRocket;
    else
        b = b + 1;
    end
    n = n + 1;
end

%% Analysis
% Initialize arrays for parameters of interest
W_wet = zeros(g,1);
W_dry = zeros(g,1);
W_prop = zeros(g,1);

L = zeros(g,1);
L_nose = zeros(g,1);
L_body = zeros(g,1);
D = zeros(g,1);

mat_body = zeros(g,1);
mat_ox_t = zeros(g,1);
mat_CC = zeros(g,1);

z_max = zeros(g,1);
ORS = zeros(g,1);
MOS = zeros(g,1);

% Extract data
for i = 1:g
    W_wet(i) = results(i).rocket.weight.total.W_wet;
    W_dry(i) = results(i).rocket.weight.total.W_dry;    
    W_prop(i) = results(i).rocket.weight.total.W_propellant;
    
    L(i) = results(i).rocket.geo.total.L;
    L_nose(i) = results(i).rocket.geo.nose.L;
    L_body(i) = results(i).rocket.geo.body.L;
    D(i) = results(i).rocket.geo.body.D;
    
    mat_body(i) = results(i).rocket.geo.body.material.no;
    mat_ox_t(i) = results(i).rocket.geo.ox_t.material.no;
    mat_CC(i) = results(i).rocket.geo.CC.material.no;
    
    z_max(i) = results(i).rocket.data.performance.z_max;
    ORS(i) = results(i).rocket.data.performance.ORS;
    MOS(i) = results(i).rocket.data.performance.MOS;
end

% Sort by weight
[W_wet, sort_index] = sort(W_wet);

W_dry = W_dry(sort_index);
W_prop = W_prop(sort_index);
L = L(sort_index);
L_nose = L_nose(sort_index);
L_body = L_body(sort_index);
D = D(sort_index);
mat_body = mat_body(sort_index);
mat_ox_t = mat_ox_t(sort_index);
mat_CC = mat_CC(sort_index);
z_max = z_max(sort_index);
ORS = ORS(sort_index);
MOS = MOS(sort_index);

% Plot
toc