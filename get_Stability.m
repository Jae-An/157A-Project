function [rocket] = get_Stability(rocket)
% Finds fin sizing for minimum 1.25 stability margin. This method assumes a
% clipped delta planform with C = C_R = 2*C_T = SS. The fin section is
% square with max thickness = 0.06C. Also finds MOS, recalculates CG using
% this fin.

body = rocket.geo.body;
fins = rocket.geo.fins;
weight = rocket.weight;

%% Solver Parameters
CNa_n = 2;
CNa_b = 4*body.L*deg2rad(12.2)/(pi*body.D); % 12.2 deg is max ORAOA
X_n = 0.666*rocket.geo.nose.L; % assume cone
X_b = body.x + body.L/2;
sumCNa_nb = CNa_n + CNa_b;
sumCNaX_nb = CNa_n*X_n + CNa_b*X_b;

X_CP = rocket.weight.total.CG_wet + rocket.data.requirements.MOS*body.D;
L = rocket.geo.total.L;
if fins.N == 3
    beta = 13.85;
elseif fins.N == 4
    beta = 16;
else
    error('Fin number wrong')
end
B = beta / ((1 + sqrt(26)/3)*body.D^2);

%% Solve for C
% p = [7*B, 12*B*(X_CP - L), 0, 12*(X_CP*sumCNa_nb - sumCNaX_nb)]; % Roots of this polynomial give possible root chords
% p_roots = roots(p); % Solve polynomial
% p_roots = p_roots(p_roots>0); % Filter out negative results
% C = min(p_roots);

%% Fill fin parameters
fins.RC = 7.5/12;%C;
fins.TC = 3.75/12;%0.5*C;
fins.SS = 4.125/12;%C;
fins.SL = 7.25/12;%0.5*C;
fins.MC = sqrt(fins.SS^2 + (fins.SL+fins.TC/2 - fins.RC/2)^2);
fins.t = 0.06*fins.RC+fins.TC/2;
fins.x = L - fins.RC; % fins located at end of rocket

a_coef = (fins.TC - fins.RC)/fins.SS;
b_coef = fins.RC;
V_fins = 4 * 5.685 / 12^3;%4 * 0.0042 * ((1/3)*(a_coef^2)*(fins.SS^3) + b_coef*fins.SS^2 + b_coef^2*fins.SS);
weight.fins.W = V_fins * 111.24; % [lb] Uses carbon fiber density
weight.fins.CG = fins.x + (11/18 * fins.RC);

%% Recalculate CG, CP, MOS
weight.total.W_wet = weight.nose.W + weight.body.W + weight.fins.W...
                   + weight.ox_t.W + weight.fuel.W + weight.ox.W + sum(weight.misc.W);
weight.total.W_dry = weight.total.W_wet - weight.total.W_propellant; 
    
W_arr = [weight.nose.W, weight.body.W, weight.fins.W, weight.payload.W, weight.ox_t.W, weight.CC.W, sum(weight.misc.W)];    
CG_arr = [weight.nose.CG, weight.body.CG, weight.fins.CG, weight.payload.CG, weight.ox_t.CG, weight.CC.CG, sum(weight.misc.W.*weight.misc.CG)/sum(weight.misc.W)]; 
weight.total.CG_dry = sum(W_arr.*CG_arr) / sum(W_arr);
weight.total.CG_wet = (weight.total.CG_dry*weight.total.W_dry + weight.fuel.CG*weight.fuel.W + weight.ox.CG*weight.ox.W) / weight.total.W_wet;

CNa_f = (1 + (body.D/2)/(fins.SS+(body.D/2))) * (4*fins.N*(fins.SS/body.D)^2 / (1+sqrt(1+(2*fins.MC/(fins.RC+fins.TC))^2)));
X_f = fins.x ...
    + fins.SL * (fins.RC + 2*fins.TC)/(3*(fins.RC+fins.TC)) ...
    + (fins.RC + fins.TC - fins.RC*fins.TC/(fins.RC+fins.TC))/6; 

X_CP = (sumCNaX_nb + CNa_f*X_f) / (sumCNa_nb + CNa_f);



%% Send data back
rocket.data.performance.MOS_wet = (X_CP - weight.total.CG_wet) / body.D;
rocket.data.performance.MOS_dry = (X_CP - weight.total.CG_dry) / body.D;
rocket.geo.body = body;
rocket.geo.fins = fins;
rocket.weight = weight;
end