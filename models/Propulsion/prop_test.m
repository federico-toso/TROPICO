% setup
clear all
close all
clc

% addpath(genpath(cd));
datapoints = 30;
minalt = 0;
maxalt = 60e3;
minmach = 0;
maxmach = 10;
throttle = 1;
model1=@prop_simpleAB;
model2=@prop_RELrocket;

scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)/2-10 50 scrsz(3)/2-10 scrsz(4)-135],'Color',[1 1 1])

th1 = zeros (datapoints,datapoints);     th2 = th1;      isp1=th1;       isp2=isp1;
altitude = linspace(minalt,maxalt,datapoints);
Mach = linspace(minmach,maxmach,datapoints);

for i = 1:datapoints
    for j = 1:datapoints
        [press, ~, ~, ~] = atmo_ISA(altitude(i));
        [th1(i,j), mp1] = model1 (throttle, Mach(j), press, 0,0);
        isp1(i,j)=th1(i,j)/(mp1*9.80665);
        [th2(i,j), mp2] = model2 (throttle, Mach(j), press, 0,0);
        isp2(i,j)=th2(i,j)/(mp2*9.80665);
    end
end

subplot(2,1,1)
hold on
colormap([1 0 0; 0 0 1])
% doublesurf(altitude/1000,Mach,isp1,hot(datapoints),altitude/1000,Mach,isp2,hot(datapoints))
surf(Mach,altitude/1000,isp1,zeros(size(isp1)))
surf(Mach,altitude/1000,isp2,ones(size(isp1)))
ylabel ('Altitude [km]', 'FontSize', 20)
xlabel ('Mach', 'FontSize', 20)
zlabel ('Isp [s]', 'FontSize', 20)
zlim(zlim+10)
view(45,30)
% caxis([0, 1])
% shading interp

subplot(2,1,2)
hold on
surf(Mach,altitude/1000,th1, zeros(size(th1)))
surf(Mach,altitude/1000,th2, ones(size(th1)))
ylabel ('Altitude [km]', 'FontSize', 20)
xlabel ('Mach', 'FontSize', 20)
zlabel ('Thrust [N]', 'FontSize', 20)
zlim(zlim+10)
view(45,30)
% caxis([0, 1])
% shading interp

% myaa(8)
disp(['min thrust,',func2str(model1),': ',num2str(min(min(th1))/1000),' kN'])
disp(['min thrust,',func2str(model2),': ',num2str(min(min(th2))/1000),' kN'])

disp(['max thrust,',func2str(model1),': ',num2str(max(max(th1))/1000),' kN'])
disp(['max thrust,',func2str(model2),': ',num2str(max(max(th2))/1000),' kN'])

disp(['max Isp,',func2str(model1),': ',num2str(max(max(isp1))),' s'])
disp(['max Isp,',func2str(model2),': ',num2str(max(max(isp2))),' s'])