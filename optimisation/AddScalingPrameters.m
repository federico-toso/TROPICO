function [LB, UB] = AddScalingPrameters (LB, UB, param)

if isfield(param,'Scale_S1_Engine')
    if param.Scale_S1_Engine(1)-param.Scale_S1_Engine(2)~=0
        LB = [LB, param.Scale_S1_Engine(1)];
        UB = [UB, param.Scale_S1_Engine(2)];
    end
end
if isfield(param,'Scale_S2_Engine')
    if param.Scale_S2_Engine(1)-param.Scale_S2_Engine(2)~=0
        LB = [LB, param.Scale_S2_Engine(1)];
        UB = [UB, param.Scale_S2_Engine(2)];
    end
end
if isfield(param,'Scale_S1_Surf')
    if param.Scale_S1_Surf(1)-param.Scale_S1_Surf(2)~=0
        LB = [LB, param.Scale_S1_Surf(1)];
        UB = [UB, param.Scale_S1_Surf(2)];
    end
end
if isfield(param,'Scale_S2_Surf')
    if param.Scale_S2_Surf(1)-param.Scale_S2_Surf(2)~=0
        LB = [LB, param.Scale_S2_Surf(1)];
        UB = [UB, param.Scale_S2_Surf(2)];
    end
end

end


