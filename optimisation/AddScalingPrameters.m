function [c0, LB, UB, phases] = AddScalingPrameters (c0, LB, UB, param, phases)

if isfield(param,'Scale_S1_Engine')
    if param.Scale_S1_Engine(1)-param.Scale_S1_Engine(2)~=0
        c0 = [c0, mean(param.Scale_S1_Engine)];
        LB = [LB, param.Scale_S1_Engine(1)];
        UB = [UB, param.Scale_S1_Engine(2)];
    end
end
if isfield(param,'Scale_S2_Engine')
    if param.Scale_S2_Engine(1)-param.Scale_S2_Engine(2)~=0
        c0 = [c0, mean(param.Scale_S2_Engine)];
        LB = [LB, param.Scale_S2_Engine(1)];
        UB = [UB, param.Scale_S2_Engine(2)];
    end
end
if isfield(param,'Scale_S1_Surf')
    if param.Scale_S1_Surf(1)-param.Scale_S1_Surf(2)~=0
        c0 = [c0, mean(param.Scale_S1_Surf)];
        LB = [LB, param.Scale_S1_Surf(1)];
        UB = [UB, param.Scale_S1_Surf(2)];
    end
end
if isfield(param,'Scale_S2_Surf')
    if param.Scale_S2_Surf(1)-param.Scale_S2_Surf(2)~=0
        c0 = [c0, mean(param.Scale_S2_Surf)];   
        LB = [LB, param.Scale_S2_Surf(1)];
        UB = [UB, param.Scale_S2_Surf(2)];
    end
end

end


