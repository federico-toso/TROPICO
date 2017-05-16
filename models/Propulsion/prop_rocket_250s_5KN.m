function [T, mp] = prop_rocket_450s_5KN(throttle, ~, P_Air,~,~)

% 0.5 m radius nozzle
Ae = pi * (1/2)^2;
% hydrazine Isp, square law losses for throttle
Isp = 450 *( -(throttle-1)^2+1);

mp = 2.04 * throttle;
T = max(0,(9.80665 * mp * Isp) - P_Air * Ae);