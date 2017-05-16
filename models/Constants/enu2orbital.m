function [a,e,i] = enu2orbital (x)

altitude   = x(1);
velocity   = x(2);
flightpath  = x(3);
heading    = x(4);
latitude   = x(5);
longitude  = x(6);
loadconstants
E = wgs84Ellipsoid;

% compute position
[X,Y,Z] = enu2ecef(0,0,0,latitude,longitude,altitude,E,'radians');
r = norm([X,Y,Z]);

% compute velocity vector
uEast   = velocity * cos(flightpath) * sin(heading) + const.wE * r * cos(latitude);
vNorth  = velocity * cos(flightpath) * cos(heading);
wUp     = velocity * sin(flightpath);
[U,V,W] = enu2ecefv(uEast,vNorth,wUp,latitude,longitude,'radians');
v = norm([U,V,W]);

% compute orbital parameters
h = cross([X,Y,Z],[U,V,W]);                     % specific angular momentum
E = (v^2)/2 - const.mu/r;                       % specific energy
a = -const.mu/(2*E);                                                        % semimayor axis
e = sqrt(abs(1-norm(h)^2/(a*const.mu)));                                    % eccentricity
i = acos(h(3)/norm(h));                                                     % inclination
% raan = atan2(h(1),-h(2));                                                   % right ascension of ascending node
% ni = acos((a*(1-e^2)-r)/(e*r));                 % true anomaly
% arg_peri = atan2(Z/sin(i),(X*cos(raan)+Y*sin(raan))) - ni;                  % argument of periapsis
% EA = 2 * atan2(sqrt((1-e)/(1+e))*tan(ni/2));    % Eccentric Anomaly
% n = sqrt(const.mu/a^3);
% T = t-1/n*(EA-e*sin(EA));                                                   % time of periapsis passage