

%% Constants for Earth ellipsoidal model 
% IERS (2003)[3]	aE=6378136.6	 bE=6356751.9	  flattening=298.25642
% const.aE = 6378136.6;  % Equatorial radius (m) of Earth
% const.bE = 6356751.9;  % Polar radius (m) of Earth
% const.eE = 0.0818192206096493;  % sqrt(1-(const.bE/const.aE)^2)

% Enhanced WGS-84 (1984)	6,378,137	6,356,752.3142	298.257223563	Global GPS
const.aE = 6378137.00;  % Equatorial radius (m) of Earth
const.bE = 6356752.3142;  % Polar radius (m) of Earth
const.eE = 0.0818191909289069;  %sqrt(1-(const.bE/const.aE)^2)  Eccentricity of Earth ellipsoid
const.flatE = 1/298.257223563;  %Earth flattening
const.wE = 7.292115e-5;  % Angular rotation of the Earth (rad/s)
% const.rE = const.aE/sqrt(1-const.eE^2*sin(pi/8)^2); %at nominal latitude 22.5

const.rE = 6375253;             % Radius of the Earth (m)
% const.wE = 7.292118e-05;        % Angular rotation of the Earth (rad/s)

%% Gravitational constants
const.g0 = 9.80665;             % Gravitational acceleration at sea level (m/s2)
const.mu = 3.986004418e14;  % Geocentric constant of gravitation (GM) [m3/s2], ref: IERS Numerical Standards

const.J2 = 1.0826359e-3;  %Second degree term in Earth's gravity potential, ref: IERS Numerical Standards
const.J3=2.532153e-7;
const.J4=1.6109876e-7;

%% Astronomic constants
const.AU = 149597870.691e3;  % Astronomical Unit (AU) [m]

%% Physical constants
const.c = 299792458; % Speed of light in the vacuum [m/s]
const.Tabs = 273.15; % Absolute temperature [K], ref: IUPAC
const.Rgas = 8.3144621;  % Universal gas constant [J/(K.mol)] ref: IUPAC
const.kB = 1.3806488e-23;  % Boltzmann constant [J/K], ref: CODATA
const.sb = 5.670373e-8; % Stefan-Boltzmann constant [W/m2.K4], ref: CODATA

%% Propulsion constants
%load('surrogate_hyperion.mat')
%const.surrogate_ejector = surrogate_ejector;
%const.surrogate_ramjet = surrogate_ramjet;
%const.surrogate_rocket = surrogate_rocket;

