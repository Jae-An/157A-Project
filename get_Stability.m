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
p = [7*B, 12*B*(X_CP - L), 0, 12*(X_CP*sumCNa_nb - sumCNaX_nb)]; % Roots of this polynomial give possible root chords
p_roots = roots(p); % Solve polynomial
p_roots = p_roots(p_roots>0); % Filter out negative results
C = min(p_roots);

%% Fill fin parameters
fins.RC = C;
fins.TC = 0.5*C;
fins.SS = C;
fins.SL = 0.5*C;
fins.MC = sqrt(1.0625)*C;
fins.t = 0.06*C;
fins.x = L - C; % fins located at end of rocket

V_fins = fins.N*(0.5*(fins.RC+fins.TC)*fins.SS)*fins.t;
weight.fins.W = V_fins * 111.24; % [lb] Uses carbon fiber density
weight.fins.CG = fins.x + (11/18 * C);

%% Recalculate CG, CP, MOS
weight.total.W_wet = weight.nose.W + weight.body.W + weight.fins.W...
                   + weight.ox_t.W + weight.fuel.W + weight.ox.W + sum(weight.misc.W);
weight.total.W_dry = weight.total.W_wet - weight.total.W_propellant; 
    
W_arr = [weight.nose.W, weight.body.W, weight.fins.W, weight.payload.W, weight.ox_t.W, weight.CC.W, weight.fuel.W, sum(weight.misc.W)];    
CG_arr = [weight.nose.CG, weight.body.CG, weight.fins.CG, weight.payload.CG, weight.ox_t.CG, weight.CC.CG, weight.fuel.CG, sum(weight.misc.W.*weight.misc.CG)/sum(weight.misc.W)]; 
weight.total.CG = sum(W_arr.*CG_arr) / sum(W_arr);

CNa_f = (1 + (body.D/2)/(fins.SS+(body.D/2))) * (4*fins.N*(fins.SS/body.D)^2 / (1+sqrt(1+(2*fins.MC/(fins.RC+fins.TC))^2)));
X_f = fins.x ...
    + fins.SL * (fins.RC + 2*fins.TC)/(3*(fins.RC+fins.TC)) ...
    + (fins.RC + fins.TC - fins.RC*fins.TC/(fins.RC+fins.TC))/6; 

X_CP = (sumCNaX_nb + CNa_f*X_f) / (sumCNa_nb + CNa_f);



%% Send data back
rocket.data.performance.MOS = (X_CP - weight.total.CG) / body.D;
rocket.geo.body = body;
rocket.geo.fins = fins;
rocket.weight = weight;
end