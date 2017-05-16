% noflyzone visualization
close

param.no_fly = [43.1*pi/180,  17.1*pi/180, 2.30e6;...  % europe
                65*pi/180,   -19.2*pi/180, 0.27e6;...  % iceland
                64.2*pi/180,  27.2*pi/180, 0.80e6;...  % norway
                18.5*pi/180,  -0.2*pi/180, 1.88e6;...  % africa
                ];
            
const.rE = 6.371e6;
            
            
%   3D plot    spherical earth
figure
[xnfz,ynfz,znfz]=sphere(100);
Xt=const.rE*xnfz;
Yt=const.rE*ynfz;
Zt=const.rE*znfz;
hold on
axis equal
rotate3d
% oldFolder = cd('other files');
I=imread('earth_blb.png');
% cd(oldFolder);
warp(Xt,Yt,Zt,I)
view(95,30)
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
% draw no fly circles --------------------------------------------------
no_fly_arc = param.no_fly(:,3)./const.rE;
[nfzlat,nfzlon] = scircle1(param.no_fly(:,1),param.no_fly(:,2),no_fly_arc);
theta=-pi:pi/100:pi;
xcirc = cos(theta);
ycirc = sin(theta);
zcirc = 0*theta;
for ici = 1 : size(nfzlat,2)
    [xnfz,ynfz,znfz] =  enu2ecef(param.no_fly(ici,3).*xcirc,param.no_fly(ici,3).*ycirc,zcirc,...
                        param.no_fly(ici,1),param.no_fly(ici,2),const.rE*cos(no_fly_arc(ici)),...
                        referenceSphere,'radians');
    plot3(xnfz,-ynfz,znfz, 'color', 'r', 'LineWidth', 1)
end