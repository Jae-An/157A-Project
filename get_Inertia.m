function [rocket] = get_Inertia(rocket)
    % ACCOUNT FOR NOZZLE
    % Setup
    W = rocket.weight;
    geo = rocket.geo;

    % Nose
    W.nose.I_self = W.nose.W/4 * ((geo.nose.D/2)^2 + 2*geo.nose.L^2);
    W.nose.I_wet = (W.nose.W*(W.total.CG_wet-W.nose.CG)^2) + W.nose.I_self;
    W.nose.I_dry = (W.nose.W*(W.total.CG_dry-W.nose.CG)^2) + W.nose.I_self;
    
    % Body
    W.body.I_self = (W.body.W/12) * (3*((geo.body.D/2)^2 + (geo.body.D/2-geo.body.t)^2) + geo.body.L^2);
    W.body.I_wet = (W.body.W*(W.total.CG_wet-W.body.CG)^2) + W.body.I_self;
    W.body.I_dry = (W.body.W*(W.total.CG_dry-W.body.CG)^2) + W.body.I_self;    
    
    % Tail (empty for now)
    
    % Fins (pt mass approx.)
    W.fins.I_self = 0;
    W.fins.I_wet = (W.fins.W*(W.total.CG_wet-W.fins.CG)^2) + W.fins.I_self;
    W.fins.I_dry = (W.fins.W*(W.total.CG_dry-W.fins.CG)^2) + W.fins.I_self;
    
    % Payload (cylinder)
    W.payload.I_self = (W.payload.W/12) * (3*(geo.payload.D/2)^2 + geo.payload.L^2);
    W.payload.I_wet = (W.payload.W*(W.total.CG_wet-W.payload.CG)^2) + W.payload.I_self;
    W.payload.I_dry = (W.payload.W*(W.total.CG_dry-W.payload.CG)^2) + W.payload.I_self;
    
    % Misc (pt mass approx.)
    W.misc.I_self = 0;
    W.misc.I_wet = (W.misc.W.*(W.total.CG_wet-W.misc.CG).^2) + W.misc.I_self;
    W.misc.I_dry = (W.misc.W.*(W.total.CG_dry-W.misc.CG).^2) + W.misc.I_self;
    
    % Propulsion
    W.ox_t.I_self = (W.ox_t.W/12) * (3*((geo.ox_t.D/2)^2 + (geo.ox_t.D/2-geo.ox_t.t)^2) + geo.ox_t.L^2);
    W.ox_t.I_wet = (W.ox_t.W.*(W.total.CG_wet-W.ox_t.CG).^2) + W.ox_t.I_self;
    W.ox_t.I_dry = (W.ox_t.W.*(W.total.CG_dry-W.ox_t.CG).^2) + W.ox_t.I_self;

    W.CC.I_self = (W.CC.W/12) * (3*((geo.CC.D/2)^2 + (geo.CC.D/2-geo.CC.t)^2) + geo.CC.L^2);
    W.CC.I_wet = (W.CC.W.*(W.total.CG_wet-W.CC.CG).^2) + W.CC.I_self;
    W.CC.I_dry = (W.CC.W.*(W.total.CG_dry-W.CC.CG).^2) + W.CC.I_self;

    W.ox.I_self = (W.ox.W/12) * (3*(geo.ox.D/2)^2 + geo.ox_t.L^2);
    W.ox.I = (W.ox.W.*(W.total.CG_wet-W.ox.CG).^2) + W.ox.I_self;

    W.fuel.I_self = (W.fuel.W/12) * (3*((geo.fuel.D_o/2)^2 + (geo.fuel.D_i/2)^2) + geo.ox_t.L^2);
    W.fuel.I = (W.fuel.W.*(W.total.CG_wet-W.fuel.CG).^2) + W.ox.I_self;

    % Totals
    W.total.I_wet = W.nose.I_wet + W.body.I_wet + W.fins.I_wet ...
                  + W.payload.I_wet + W.ox_t.I_wet + W.CC.I_wet + sum(W.misc.I_wet) ...
                  + W.ox.I + W.fuel.I;
    W.total.I_dry = W.nose.I_dry + W.body.I_dry + W.fins.I_dry ...
                  + W.payload.I_dry + W.ox_t.I_dry + W.CC.I_dry + sum(W.misc.I_dry);
              
    %% Return
   rocket.geo = geo;
   rocket.weight = W;
end