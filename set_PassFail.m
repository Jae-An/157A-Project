function [is_Good] = set_PassFail(rocket)
    
    % Need to add geo check (see if internal components fit in length of aeroshell)
    geo = rocket.geo;
    dx = 2/12;
    L_internal = geo.payload.x + geo.payload.L + sum(geo.misc.L) + 3*dx ...
               + geo.ox_t.L + geo.CC.L;
    
    if rocket.data.performance.ORS <= rocket.data.requirements.ORS
        % Off-rail speed check
        is_Good = false;
   
    elseif rocket.data.performance.z_max < rocket.data.requirements.z_max
        % Apogee check
        is_Good = false;
    
    elseif rocket.data.performance.MOS < rocket.data.requirements.MOS
        % Margin of stability check
        is_Good = false;
        
    elseif (geo.total.L < L_internal) ...
        || ((geo.total.L-L_internal)/geo.total.L < 0.05)
        % Geometry check (internals do not fit or there isn't enough
        % margin)
        is_Good = false;
        
    else
        is_Good = true;

    end

end