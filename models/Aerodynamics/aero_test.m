% setup
% close all
% clear all
% clc
%{
% addpath(genpath('C:\Users\ypb14182\Desktop\IACrob'));
% datapoints = 50;
% altitude = 1e3;
% minalp = deg2rad(-10);
% maxalp = deg2rad(50);
% minmach = 0.5;
% maxmach = 1;
% model1=@aero_FSPLUK;
% vehicle = load_vehicle('skylon');
% 
% CL = zeros (datapoints,datapoints);     CD = CL;
% Mach = linspace(minmach,maxmach,datapoints);
% alpha = linspace(minalp,maxalp,datapoints);
% 
% [~, dens, ss, ~] = atmo_ISA(altitude);
% 
% for i = 1:datapoints
%     for j = 1:datapoints
% %         [CD(i,j), CL(i,j)] = model1 (Mach(i), deg2rad(alpha(j)), vehicle, 2);
%         [CD(i,j), CL(i,j)] = model1 (Mach(i), alpha(j), vehicle, 2, 2);
%     end
% end
% 
% [Mgrid, agrid] = ndgrid (Mach(:), alpha(:));
% CLI = griddedInterpolant(Mgrid, agrid, CL,'cubic');
% CDI = griddedInterpolant(Mgrid, agrid, CD);
% 
% refine = 3;
% MachD=linspace(minmach-1,maxmach+1,datapoints*refine);
% alphaD = linspace(minalp-0.2,maxalp+0.2,datapoints*refine);
% for i = 1:datapoints*refine
%     for j = 1:datapoints*refine
%         CLID(i,j)=CLI(MachD(i), alphaD(j));
%         CDID(i,j)=CDI(MachD(i), alphaD(j));
%     end
% end
% 
% [MDP,aDP]=meshgrid(MachD,alphaD);
% 
% scrsz = get(groot,'ScreenSize');
% figure('Position',[scrsz(3)/2-10 50 scrsz(3)/2-10 scrsz(4)-135],'Color',[1 1 1])
% 
% subplot(2,1,1)
% hold on
% surf(alpha,Mach,CL,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% mesh(alphaD,MachD,CLID)
% xlabel ('alpha [deg]', 'FontSize', 20)
% ylabel ('Mach', 'FontSize', 20)
% zlabel ('CL', 'FontSize', 20)
% % zlim(zlim+10)
% view(45,30)
% % shading interp
% 
% subplot(2,1,2)
% hold on
% surf(alpha,Mach,CD,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% mesh(alphaD,MachD,CDID)
% xlabel ('alpha [deg]', 'FontSize', 20)
% ylabel ('Mach', 'FontSize', 20)
% zlabel ('CD', 'FontSize', 20)
% % zlim(zlim+10)
% view(45,30)
% % shading interp
% 
% % myaa(8)
% 
% figure, hold on
% surf(alpha,Mach,CL./CD,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% xlabel ('alpha [deg]', 'FontSize', 20)
% ylabel ('Mach', 'FontSize', 20)
% zlabel ('L/D', 'FontSize', 20)
% view(45,30)
%}
  %%  
datapoints = 1000;
altitude = 1e3;
minalp = deg2rad(-10);
maxalp = deg2rad(30);
minmach = 0;
maxmach = 6.5;
model1=@aero_X34;
% vehicle = load_vehicle('AIAAHYP17');

CL = zeros (datapoints,datapoints);     CD = CL;
Mach = linspace(minmach,maxmach,datapoints);
alpha = linspace(minalp,maxalp,datapoints);

[~, dens, ss, ~] = atmo_ISA(altitude);
%specific for x34
CDdata = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',2 ,'A3:C149');
CLdata = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',1 ,'A4:C98');
CDdata(~any(~isnan(CDdata), 2),:)=[];
CLdata(~any(~isnan(CLdata), 2),:)=[];

for i = 1:datapoints
    for j = 1:datapoints
        [CD(i,j), CL(i,j)] = model1 (Mach(i), alpha(j), [], CDdata, CLdata);
        alphadata(i,j) = alpha(j);
        Machdata(i,j) = Mach(i);
    end
end

% CDsub = CD(Machdata<=1);
% CLsub = CL(Machdata<=1);
% Machdatasub = Machdata(Machdata<=1);
% alphadatasub = alphadata(Machdata<=1);
% CDsup = CD(Machdata>1);
% CLsup = CL(Machdata>1);
% Machdatasup = Machdata(Machdata>1);
% alphadatasup = alphadata(Machdata>1);

scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)/2-10 50 scrsz(3)/2-10 scrsz(4)-135],'Color',[1 1 1])

subplot(2,1,1)
hold on
surf(alpha,Mach,CL,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel ('alpha [rad]', 'FontSize', 20)
ylabel ('Mach', 'FontSize', 20)
zlabel ('CL', 'FontSize', 20)
view(45,45)

subplot(2,1,2)
hold on
surf(alpha,Mach,CD,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel ('alpha [rad]', 'FontSize', 20)
ylabel ('Mach', 'FontSize', 20)
zlabel ('CD', 'FontSize', 20)
view(45,45)

figure, hold on
surf(alpha,Mach,CL./CD,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel ('alpha [rad]', 'FontSize', 20)
ylabel ('Mach', 'FontSize', 20)
zlabel ('L/D', 'FontSize', 20)
view(45,45)

