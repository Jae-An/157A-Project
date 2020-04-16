function [is_Good] = set_PassFail(rocket)
    
    if (rocket.data.performance.ORS > rocket.data.requirements.ORS)...
    && (rocket.data.performance.z_max >= rocket.data.requirements.z_max)
        
        is_Good = true;
        
    else
        
        is_Good = false;

    end

end