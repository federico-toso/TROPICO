function [T, mp] = prop_rocket_450s_1MN (throttle, ~, P_Air,~,~)

% 1 m radius nozzle
Ae = pi;
% hydrogen Isp, square law losses for throttle
Isp = 450; %%%%%%%%%%%%%%%%%%%%%%%%%%%*( -(throttle-1)^2+1);

mp = 226.60 * throttle;
T = max(0,(9.80665 * mp * Isp) - P_Air * Ae);