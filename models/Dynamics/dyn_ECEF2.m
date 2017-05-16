function [dx, con] = dyn_ECEF2(t, x, control, const, phases)
% function [dx, con] = dyn_ECEF2(t, x, control, vehicle, param, phases, const, iphase)
%{
[dx, con] = dyn_ECEF2(t, x, control, vehicle, param, phases, const, iphase)

System model including time-derivates of the equations of motion for a
point with variable mass moving along a trajectory above a spherical rotating Earth.
 
INPUT
t = time
x = [r, v, gamma, chi, lat, lon, m] = state vector containing: radius r, relative velocity v, flight path
angle gamma, flight heading angle chi, latitude lat, longitude lon and vehicle mass m
control = [time alpha throttle bank] = control vector for the net thrust vector acting on the centre of mass of the vehicle
    alpha, bank = angle of attack and bank angle, rad
    throttle = fraction of max thrust applied, [0, 1]
vehicle = vehicle-specific values such as dimensions, mass, configurations, aerodyn coeffs, etc. 
param = common mission parameters
phases = phase structure containing phase number, models, etc.
const = physical const such as the radius of the earth
iphase = phase number (integer >= 1)

OUTPUT
dx = [dr/dt, dv/dt, dgamma/dt, dchi/dt, dlat/dt, dlon/dt, dm/dt] = values of the time derivates of the state vector
con = [accx, accz] = instantaneous accelerations, used for constraints in optimisation

FUNCTIONS CALLED
1- Interpolation law, e.g. pchip
2- Atmospheric model, e.g. atmo_ISA
3- Aerodynamic model, e.g. aero_SSTO_mex
4- Propulsion model, e.g. prop_RELrocket

CHANGE LOG
- C Maddock, 03/2016, corrected EoM for chi, introduced gravity model based
on spherical harmonics

(c) 2015 Centre for Future Air-Space Transportation Technology, Univeristy of Strathclyde
%}

% State vector components
h = x(1); % altitude from the surface of the Earth to the vehicle, m 
v = x(2); % vehicle velocity relative to the Earth, m/s
fpa = x(3); % flight path angle (angle between the relative velocity and the local horizontal axis), rad
chi = x(4); % flight heading angle (angle between the velocity vector projected onto the local horizontal plane, and North), rad
lat = x(5);  % latitude (+ North), rad
lon = x(6); % longitude (+ East), rad
m = x(7); % vehicle mass, kg
con.state=x;

rE = const.aE*(1-const.flatE*sin(lat)^2);  % WSG-84 model values for spheroid
r = h + rE;
[gr, gt]=dyn_gravity(r, lat,const);

%% Control law (h-alpha)
% The control law is up to a 4 column matrix that gives time (s), angle of attack
% alpha (rad), throttle and bank angle beta (rad).
% The following interpolates this matrix to determine the value of alpha
% at the given time instant during the integration
controls = phases.fg;
[controls] = phases.inte(control,t);
controls=max(min(controls(:),phases.cbounds(:,2)), phases.cbounds(:,1));    % limit interpolation out of bounds
alpha = controls(1);
throttle = controls(2);
bank = controls(3);
con.controls = controls;

if h <= 0 || m <= 0 || v <= 0
    dx =zeros(7,1);
    con.acc=[0, 0];
    con.Mn = 0;
    con.forces = zeros(1,4);
    con.temp = 0;
    return
end

%% Atmospheric model
% to determine pressure P, density rho and the speed of sound a at a given altitude h
[P_h, rho_h, a_h, T_h] = phases.atmo(h);
Mn = abs(v)/a_h;          % mach number

%% Propulsion model
% Determines the level of thrust (N) and fuel mass flow rate (kg/s) for the
% vehicle at a given altitude and mach number
alpha_eng=alpha+phases.vehicle.Tangle;  % accounts for any thrust angular offset from body axis
[FT, mp] = phases.prop (throttle, Mn, P_h, T_h, t);

%% Aerodynamics
% Determines the coefficients of lift and drag given the mach number and
% angle of attack, and calculates the lift and drag coefficients acting on
% the vehicle
[CD, CL] =  phases.aero(Mn, alpha, phases.vehicle, phases, 1);
q = 0.5*rho_h*(a_h*Mn)^2;
L = CL*phases.vehicle.Sgross*q;
D = CD*phases.vehicle.Sgross*q;

%% Thermodynamics
[con.temp] = phases.thermal(T_h, rho_h, v, phases.vehicle, const.sb);

%% forces scaling 
FT=phases.scaling.engine*FT; 
mp=phases.scaling.engine*mp; 
L=phases.scaling.surf*L;
D=phases.scaling.surf*D;

%% Equations of motion
% Vehicle dynamics and kinematics as a function of time
dh = v.*sin(fpa);
dlat = v.*cos(fpa).*cos(chi)./r;
if  abs(lat-pi/2) < eps(1)
    dlon = 0;
else
    dlon = v.*cos(fpa).*sin(chi)./(r.*cos(lat));
end

Fx = (FT.*cos(alpha_eng)*cos(bank) - D)./m - gr.*sin(fpa) + gt*cos(fpa)*cos(chi);
dv = Fx + const.wE.^2.*r.*cos(lat).*(sin(fpa).*cos(lat) - cos(fpa).*cos(chi).*sin(lat));
if abs(v) < eps(1)
    dfpa = 0;
else
    Fz = (FT*sin(alpha_eng)+L)*cos(bank)/m - gr*cos(fpa) - gt*sin(fpa)*cos(chi);
    dfpa = (v/r)*cos(fpa) + Fz/v + (const.wE^2*r/v)*cos(lat)*(sin(fpa)*cos(chi)*sin(lat) + cos(fpa)*cos(lat)) + 2*const.wE*sin(chi)*cos(lat);
end
if (abs(lat - pi/2) < eps(1) || abs(lat + pi/2) < eps(1) || abs(fpa - pi/2) < eps(1) || abs(v) < eps(1))
    dchi = 0;
else
    Fy = (FT*sin(alpha_eng)+L)*sin(bank)/m - gt*sin(chi);
    dchi = (v/r)*cos(fpa)*sin(chi)*tan(lat) + Fy/(v*cos(fpa)) ...
        + const.wE^2*r*(sin(chi)*sin(lat)*cos(lat))/(v*cos(fpa)) + 2*const.wE*(sin(lat)-tan(fpa)*cos(chi)*cos(lat));
end

dm = -mp;

dx = [dh; dv; dfpa; dchi; dlat; dlon; dm];
% dx,x,keyboard

accx = (FT.*cos(alpha_eng)-D.*cos(alpha)+L.*sin(alpha))./m;%-gr*sin(fpa+alpha);
accz = -(L.*cos(alpha)+D.*sin(alpha)+FT.*sin(alpha_eng))./m;%+gr*cos(fpa+alpha);

con.acc = [accx, accz];
% con.state = x;
% con.control = controls;
con.Mn = Mn;
con.forces = [L,D,FT,q];
return

function [gr, gt]=dyn_gravity(r, lat,const)
% loadconstants

% g  = const.g0*(const.rE/r)^2;  % value of gravity at a given altitude
% gr = g; gt = 0;

phi=pi/2-lat;
gr=(const.mu/r^2*(1 - 1.5*const.J2*(3*cos(phi)^2-1)*(const.rE/r)^2 ...
    -2*const.J3*cos(phi)*(5*cos(phi)^2-3)*(const.rE/r)^3 ...
    -(5/8)*const.J4*(35*cos(phi)^4-30*cos(phi)^2+3)*(const.rE/r)^4));
gt=-3*const.mu*sin(phi)*cos(phi)*(const.rE/r)^2*(const.J2+0.5*const.J3*(5*cos(phi)^2-1)...
    *(const.rE/r)/cos(phi)+(5/6)*const.J4*(7*cos(phi)^2-1)*(const.rE/r)^2)/r^2;
