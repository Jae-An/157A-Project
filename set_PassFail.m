function [is_Good] = set_PassFail(rocket)
    
    % Need to add geo check (see if internal components fit in length of aeroshell)
    geo = rocket.geo;
    dx = 1/12;
    L_internal = geo.payload.x + geo.payload.L + sum(geo.misc.L) + 3*dx ...
               + geo.ox_t.L + geo.CC.L;
    
    if rocket.data.performance.v_ORS <= rocket.data.requirements.v_ORS
        % Off-rail speed check
        is_Good = false;
   
    elseif rocket.data.performance.z_max < rocket.data.requirements.z_max
        % Apogee check
        is_Good = false;
    
    elseif rocket.data.performance.MOS_dry < rocket.data.requirements.MOS ...
        || rocket.data.performance.MOS_wet < rocket.data.requirements.MOS ...
        || ~isreal([rocket.data.performance.MOS_dry, rocket.data.performance.MOS_wet])
        % Margin of stability check
        is_Good = false;
        
    elseif (geo.total.L < L_internal) || ((geo.total.L-L_internal) < 0.5)
        % Geometry check (internals do not fit or there isn't enough
        % margin)
        is_Good = false;
        
    else
        is_Good = true;

    end

end