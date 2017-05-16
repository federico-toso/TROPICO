function [T, mp] = prop_IAC16_TUG (throttle, ~, P_Air,~,~)

Ae = pi * (1/2)^2;
Isp= 250;
g0=9.81;
Tvac=500*g0;
T  = (Tvac - P_Air*Ae)*throttle;
mp = (Tvac/(g0*Isp))*throttle;