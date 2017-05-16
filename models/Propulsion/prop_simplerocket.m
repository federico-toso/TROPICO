function [T, mp] = prop_simplerocket (throttle, ~, P_Air,~,~)

%{
rocket engine roughly modeled on REL's SABRE data

INPUT
-throttle = value from 0 to 1
-P_Air = air pressure [Pa]
-Mach = Mach number (unused)

OUTPUT
-T = thrust, [N]
-mp = fuel mass flow [kg/s]

(c) 2015, F Toso, Centre for Future Air-Space Transportation
Technology, Univeristy of Strathclyde
%}


AreaExit    = 3.4  ;
TsingleRoc  = 1e6 * throttle;
SpecImpRoc  = 450*9.81 ;     
MixRatioRoc = 6.0    ;

MdotFuel = ( 4 * TsingleRoc / SpecImpRoc )/( 1 + MixRatioRoc );
MdotOxid = MdotFuel * MixRatioRoc;

T = 4 * (TsingleRoc - AreaExit * P_Air);
T = max(min(T,1e8),0);
mp=MdotFuel+MdotOxid;