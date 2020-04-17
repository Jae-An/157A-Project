function [is_Good] = set_PassFail(rocket)
    
    % Need to add geo check (see if internal components fit in length of aeroshell)
    
    if (rocket.data.performance.ORS > rocket.data.requirements.ORS)...
    && (rocket.data.performance.z_max >= rocket.data.requirements.z_max)...
    && (rocket.data.performance.MOS >= rocket.data.requirements.MOS)   
        is_Good = true;
        
    else
        
        is_Good = false;

    end

end