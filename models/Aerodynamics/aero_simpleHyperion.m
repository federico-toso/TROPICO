function [CD, CL] =  aero_simpleHyperion(~, alpha, ~, ~, ~)

%{
aero_simpleHyperion: Function to obtain aerodynamic coefficients

Inputs:
* alpha : angle of attack [deg]

Outputs:
* CD : Drag coefficient
* CL : Lift coefficient

Author: Federico Toso
email: federico.toso@strath.ac.uk
%}

Cla=0.012;
k=2.09;
Cd0=0.012;
CL=Cla*(rad2deg(alpha)+1.14);
CD=Cd0+k*CL^2;