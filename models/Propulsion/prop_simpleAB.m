function [T, mp] = prop_simpleAB (throttle, Mach, P_Air,~,~)

%{
the airbreathing engine is here  roughly modeled on REL's SABRE data

INPUT
-throttle = value from 0 to 1
-P_Air = air pressure [Pa]
-Mach = Mach number

OUTPUT
-T = thrust, [N]
-mp = fuel mass flow [kg/s]

(c) 2015, F Toso, Centre for Future Air-Space Transportation
Technology, Univeristy of Strathclyde
%}

nacelle = 2;
mp = throttle*25*nacelle;
% height penalty --------------------------------------------------
% Isp = 6500 SL, 0 vacuum
Isp = 8000*(P_Air/50e3)^0.2;
% decrease in Isp due to Mach, parabolic, peak@Mach = 3.5
Isp = Isp*(-0.043*Mach^2+0.3*Mach+0.45);
% decrease in Isp due to throttle ---------------------------------
% 100% at full throttle
Isp = Isp*(-0.2*throttle^2+0.4*throttle+0.8);
% thrust = massflow*g*Isp -----------------------------------------
T = mp*9.80665*Isp;