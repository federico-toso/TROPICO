function [phases] = ExtractScalingPrameters (param, c0, phases)

Scale_S2_Surf = 1;
Scale_S1_Surf = 1;
Scale_S2_Engine = 1;
Scale_S1_Engine = 1;

if param.ScalingActive == 1
    k=length(c0);
    if isfield(param,'Scale_S2_Surf')
        if param.Scale_S2_Surf(1)-param.Scale_S2_Surf(2)==0
            Scale_S2_Surf = param.Scale_S2_Surf(1);
        else
            Scale_S2_Surf = c0(k);
            k = k-1;
        end
    end
    if isfield(param,'Scale_S1_Surf')
        if param.Scale_S1_Surf(1)-param.Scale_S1_Surf(2)==0
            Scale_S1_Surf = param.Scale_S1_Surf(1);
        else
            Scale_S1_Surf = c0(k);
            k = k-1;
        end
    end
    if isfield(param,'Scale_S2_Engine')
        if param.Scale_S2_Engine(1)-param.Scale_S2_Engine(2)==0
            Scale_S2_Engine = param.Scale_S2_Engine(1);
        else
            Scale_S2_Engine = c0(k);
            k = k-1;
        end
    end
    if isfield(param,'Scale_S1_Engine')
        if param.Scale_S1_Engine(1)-param.Scale_S1_Engine(2)==0
            Scale_S1_Engine = param.Scale_S1_Engine(1);
        else
            Scale_S1_Engine = c0(k);
        end
    end    
end

% add noting where to use stage 1 and where to use stage 2 settings
for ip = 1:numel(phases)
    if isequal(phases(ip).vehicle,phases(1).vehicle)    % STAGE 1
        phases(ip).scaling.engine = Scale_S1_Engine;
        phases(ip).scaling.surf = Scale_S1_Surf;
    else                                                % STAGE 2+
        phases(ip).scaling.engine = Scale_S2_Engine;
        phases(ip).scaling.surf = Scale_S2_Surf;
    end
end

end