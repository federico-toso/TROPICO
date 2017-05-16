function [press, dens, sspeed, temp] = atmo_ISA(altitude)

%{
atmo_ISA: Function to obtain atmospheric characteristics with ISA 
atmospheric model

Inputs:
* altitude : self-explanatory 

Outputs:
* press : air pressure at the requested altitude [Pa]
* dens : air density at the requested altitude [kg/m3]
* sspeed : air sound speed at the requested altitude [m/s]
* temp : air temperature at the requested altitude [K]

Author: Federico Toso
email: federico.toso@strath.ac.uk
%}

% ATMOSPHERIC DATA ---------------------------
altitudes = [0 11000 20000 32000 47000 51000 71000 84852];
staticpres= [101325 22632.1 5474.89 868.02 110.91 66.94 3.96 0.3734];
standtemp = [288.15 216.65 216.65 228.65 270.65 270.65 214.65 186.87];
lapserate = [-0.0065 0 0.001 0.0028 0 -0.0028 -0.002 0];
        R = 8.31432;
       g0 = 9.80665;
        M = 0.0289644;
        i = 0;

if     altitude < altitudes(1)
    press = 101325;
     temp = 288.15;
elseif altitude <= altitudes(2)
    i=1;
elseif altitude <= altitudes(3)
    i=2;
elseif altitude <= altitudes(4)
    i=3;
elseif altitude <= altitudes(5)
    i=4;
elseif altitude <= altitudes(6)
    i=5;
elseif altitude <= altitudes(7)
    i=6;
elseif altitude <= altitudes(8)
    i=7;
end

if i==1||i==3||i==4||i==6||i==7
    press = staticpres(i)*(standtemp(i)/(standtemp(i)+lapserate(i)*(altitude-altitudes(i))))^(g0*M/(R*lapserate(i)));
    temp = standtemp(i)+lapserate(i)*(altitude-altitudes(i));
elseif i==2||i==5
    press = staticpres(i)*exp((-g0*M*(altitude-altitudes(i))/(R*standtemp(i))));
    temp = standtemp(i)+lapserate(i)*(altitude-altitudes(i));
elseif altitude > altitudes(8)
    press = 0;
    temp = 186.87;
end
dens=press/(287.058*temp);
sspeed=sqrt(1.4*287.058*temp);