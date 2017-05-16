function [CD, CL] =  aero_X34(Mn, alpha, ~, CDdata, CLdata)
% filter inputs
Mn = (Mn < 0.4) * 0.4 + (Mn >= 0.4) * (Mn <= 6) * Mn + (Mn > 6) * 6;
%{
% alpha = alpha*180/pi;
% alpha = (alpha < -4) * -4 + (alpha >= -4) * (alpha <= 24) * alpha + (alpha > 24) * 24;
% alpha = alpha/180*pi;
%}
alpha = (alpha < -0.0698) * (-0.0698) +...
        (alpha >= -0.0698) * (alpha <= 0.4189) * alpha +...
        (alpha > 0.4189) * 0.4189;
%{
% CDdata = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',2 ,'A3:C149');
% CDdata(~any(~isnan(CDdata), 2),:)=[];
FCD = scatteredInterpolant(CDdata(:,1),CDdata(:,2),CDdata(:,3));

CD = FCD(Mn,alpha);

% CLdata = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',1 ,'A4:C98');
% CLdata(~any(~isnan(CLdata), 2),:)=[];
FCL = scatteredInterpolant(CLdata(:,1),CLdata(:,2),CLdata(:,3));

CL = FCL(Mn,alpha);
%}

if Mn<=1
    %R-square: 0.9912
    CD = -0.01157 +0.08119 *Mn -0.01673 *alpha +0.08063 *Mn*alpha +2.811 *alpha^2;
    %R-square: 0.9954
    CL = 0.134 -0.000958 *Mn +2.869 *alpha;
else
    %R-square: 0.9926
    CD = 0.1311 -0.05283 *Mn +0.3676 *alpha +0.01508 *Mn^2 -0.2063 *Mn*alpha +...
        2.249 *alpha^2 -0.001473 *Mn^3 +0.01884 *Mn^2*alpha -0.04375 *Mn*alpha^2;
    %R-square: 0.9931
    CL = 0.5201 -0.4455 *Mn +3.234 *alpha +0.1182 *Mn^2 -0.5418 *Mn*alpha +...
        -0.01025 *Mn^3 +0.04499 *Mn^2*alpha;
end

% w/ sub-super separation (horrible graph)
% CL
%{
% % CL supersonic
% % load xlsx
% CLdatasup = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',1 ,'A4:C50');
% % purge NaN
% CLdatasup(~any(~isnan(CLdatasup), 2),:)=[];
% % separate columns
% Mn = CLdatasup(:,1);
% alpha = CLdatasup(:,2);
% CL = CLdatasup(:,3);
% % evaluate coefficients
% CLcoeff2 = [ones(size(CL,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CL;
% CLexperimental2 = CLcoeff2(1) + CLcoeff2(2).*Mn + CLcoeff2(3).*alpha + ...
%     CLcoeff2(4).*Mn.^2 + CLcoeff2(5).*alpha.^2 + CLcoeff2(6).*Mn.*alpha; 
% CLcoeff3 = [ones(size(CL,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
%     Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2] \ CL;
% CLexperimental3 = CLcoeff3(1) + CLcoeff3(2).*Mn + CLcoeff3(3).*alpha + ...
%     CLcoeff3(4).*Mn.^2 + CLcoeff3(5).*alpha.^2 + CLcoeff3(6).*Mn.*alpha + ...
%     CLcoeff3(7).*Mn.^3 + CLcoeff3(8).*alpha.^3 + CLcoeff3(9).*Mn.^2.*alpha + CLcoeff3(10).*Mn.*alpha.^2; 
% % visualize
% figure
% hold on
% scatter3(Mn, alpha, CL)
% % scatter3(Mn, alpha, CLexperimental2)
% % scatter3(Mn, alpha, CLexperimental3)
% % analyse
% R2= corrcoef(CL,CLexperimental2);
% R2(2,1) %.9936
% R3 = corrcoef(CL,CLexperimental3);
% R3(2,1) %.9995
% % visualize result
% testpoints = 50;
% alpha = linspace(-5,25,testpoints);
% Mn = linspace(1,6,testpoints);
% CLtest2 = zeros (testpoints,testpoints);
% CLtest3=CLtest2;
% for i = 1:testpoints
%     for j = 1:testpoints
%         CLtest2(i,j) = CLcoeff2(1) + CLcoeff2(2).*Mn(j) + CLcoeff2(3).*alpha(i) + ...
%             CLcoeff2(4).*Mn(j).^2 + CLcoeff2(5).*alpha(i).^2 +...
%             CLcoeff2(6).*Mn(j).*alpha(i);
%         CLtest3(i,j) = CLcoeff3(1) + CLcoeff3(2).*Mn(j) + CLcoeff3(3).*alpha(i) + ...
%             CLcoeff3(4).*Mn(j).^2 + CLcoeff3(5).*alpha(i).^2 +...
%             CLcoeff3(6).*Mn(j).*alpha(i)+ CLcoeff3(7).*Mn(j).^3 + ...
%             CLcoeff3(8).*alpha(i).^3 + CLcoeff3(9).*Mn(j).^2.*alpha(i) + CLcoeff3(10).*Mn(j).*alpha(i).^2; 
%     end
% end
% surf(Mn,alpha,CLtest2,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% % surf(Mn,alpha,CLtest3,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%}
%{
% CL subsonic
% load xlsx
CLdatasub = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',1 ,'A52:C98');
% purge NaN
CLdatasub(~any(~isnan(CLdatasub), 2),:)=[];
% separate columns
Mn = CLdatasub(:,1);
alpha = CLdatasub(:,2);
CL = CLdatasub(:,3);
% evaluate coefficients
CLcoeff2 = [ones(size(CL,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CL;
CLexperimental2 = CLcoeff2(1) + CLcoeff2(2).*Mn + CLcoeff2(3).*alpha + ...
    CLcoeff2(4).*Mn.^2 + CLcoeff2(5).*alpha.^2 + CLcoeff2(6).*Mn.*alpha; 
CLcoeff3 = [ones(size(CL,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
    Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2] \ CL;
CLexperimental3 = CLcoeff3(1) + CLcoeff3(2).*Mn + CLcoeff3(3).*alpha + ...
    CLcoeff3(4).*Mn.^2 + CLcoeff3(5).*alpha.^2 + CLcoeff3(6).*Mn.*alpha + ...
    CLcoeff3(7).*Mn.^3 + CLcoeff3(8).*alpha.^3 + CLcoeff3(9).*Mn.^2.*alpha + CLcoeff3(10).*Mn.*alpha.^2; 
% visualize
figure
hold on
scatter3(Mn, alpha, CL)
% scatter3(Mn, alpha, CLexperimental2)
% scatter3(Mn, alpha, CLexperimental3)
% analyse
R2= corrcoef(CL,CLexperimental2);
R2(2,1) %.9975
R3 = corrcoef(CL,CLexperimental3);
R3(2,1) %.9996
% visualize result
testpoints = 50;
alpha = linspace(-5,25,testpoints);
Mn = linspace(0,1,testpoints);
CLtest2 = zeros (testpoints,testpoints);
CLtest3=CLtest2;
for i = 1:testpoints
    for j = 1:testpoints
        CLtest2(i,j) = CLcoeff2(1) + CLcoeff2(2).*Mn(j) + CLcoeff2(3).*alpha(i) + ...
            CLcoeff2(4).*Mn(j).^2 + CLcoeff2(5).*alpha(i).^2 +...
            CLcoeff2(6).*Mn(j).*alpha(i);
        CLtest3(i,j) = CLcoeff3(1) + CLcoeff3(2).*Mn(j) + CLcoeff3(3).*alpha(i) + ...
            CLcoeff3(4).*Mn(j).^2 + CLcoeff3(5).*alpha(i).^2 +...
            CLcoeff3(6).*Mn(j).*alpha(i)+ CLcoeff3(7).*Mn(j).^3 + ...
            CLcoeff3(8).*alpha(i).^3 + CLcoeff3(9).*Mn(j).^2.*alpha(i) + CLcoeff3(10).*Mn(j).*alpha(i).^2; 
    end
end
surf(Mn,alpha,CLtest2,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
% surf(Mn,alpha,CLtest3,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%}
%{
if Mn > 1
    CLcoeff2 = [0.592865874015661,-0.427784602800927,0.0513145741288362,...
                0.0527829813941431,2.30019840030130e-05,-0.00415327890558871];
    CL = CLcoeff2(1) + CLcoeff2(2).*Mn + CLcoeff2(3).*alpha + ...
            CLcoeff2(4).*Mn.^2 + CLcoeff2(5).*alpha.^2 +...
            CLcoeff2(6).*Mn.*alpha;
else
    CLcoeff2 = [-0.0253014086504220,0.390525919458304,0.0634086128871136,...
        -0.218401797381753,-0.000201263966841061,-0.0147734362755196];
    CL = CLcoeff2(1) + CLcoeff2(2).*Mn + CLcoeff2(3).*alpha + ...
            CLcoeff2(4).*Mn.^2 + CLcoeff2(5).*alpha.^2 +...
            CLcoeff2(6).*Mn.*alpha;
end
%}
% CD
%{
% CD hypersonic
% load xlsx
CDdatasup = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',2 ,'A3:C48');
% purge NaN
CDdatasup(~any(~isnan(CDdatasup), 2),:)=[];
% separate columns
Mn = CDdatasup(:,1);
alpha = CDdatasup(:,2);
CD = CDdatasup(:,3);
% evaluate coefficients
CDcoeffx = [ones(size(CD,1),1), Mn, alpha, alpha.^2] \ CD;
CDexperimental2x = CDcoeffx(1) + CDcoeffx(2).*Mn + CDcoeffx(3).*alpha + ...
    CDcoeffx(4).*alpha.^2; 
CDcoeff2 = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CD;
CDexperimental2 = CDcoeff2(1) + CDcoeff2(2).*Mn + CDcoeff2(3).*alpha + ...
    CDcoeff2(4).*Mn.^2 + CDcoeff2(5).*alpha.^2 + CDcoeff2(6).*Mn.*alpha; 
CDcoeff3 = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
    Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2] \ CD;
CDexperimental3 = CDcoeff3(1) + CDcoeff3(2).*Mn + CDcoeff3(3).*alpha + ...
    CDcoeff3(4).*Mn.^2 + CDcoeff3(5).*alpha.^2 + CDcoeff3(6).*Mn.*alpha + ...
    CDcoeff3(7).*Mn.^3 + CDcoeff3(8).*alpha.^3 + CDcoeff3(9).*Mn.^2.*alpha + CDcoeff3(10).*Mn.*alpha.^2; 
% visualize
figure
hold on
scatter3(Mn, alpha, CD)
% scatter3(Mn, alpha, CDexperimental2)
% scatter3(Mn, alpha, CDexperimental3)
% analyse
R2= corrcoef(CD,CDexperimental2x);
R2(2,1) %.9879
R2= corrcoef(CD,CDexperimental2);
R2(2,1) %.9920
R3 = corrcoef(CD,CDexperimental3);
R3(2,1) %.9995
% visualize result
testpoints = 50;
alpha = linspace(-5,25,testpoints);
Mn = linspace(1,6,testpoints);
CDtest2 = zeros (testpoints,testpoints);
CDtest3=CDtest2;
CDtest2x = CDtest2;
for i = 1:testpoints
    for j = 1:testpoints
        CDtest2x(i,j) = CDcoeffx(1) + CDcoeffx(2).*Mn(j) + CDcoeffx(3).*alpha(i) + ...
            CDcoeffx(4).*alpha(i).^2;
        CDtest2(i,j) = CDcoeff2(1) + CDcoeff2(2).*Mn(j) + CDcoeff2(3).*alpha(i) + ...
            CDcoeff2(4).*Mn(j).^2 + CDcoeff2(5).*alpha(i).^2 +...
            CDcoeff2(6).*Mn(j).*alpha(i);
        CDtest3(i,j) = CDcoeff3(1) + CDcoeff3(2).*Mn(j) + CDcoeff3(3).*alpha(i) + ...
            CDcoeff3(4).*Mn(j).^2 + CDcoeff3(5).*alpha(i).^2 +...
            CDcoeff3(6).*Mn(j).*alpha(i)+ CDcoeff3(7).*Mn(j).^3 + ...
            CDcoeff3(8).*alpha(i).^3 + CDcoeff3(9).*Mn(j).^2.*alpha(i) + CDcoeff3(10).*Mn(j).*alpha(i).^2; 
    end
end
surf(Mn,alpha,CDtest2x,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%}
%{
% CD supersonic
% load xlsx
CDdatasup = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',2 ,'A27:C72');
% purge NaN
CDdatasup(~any(~isnan(CDdatasup), 2),:)=[];
% separate columns
Mn = CDdatasup(:,1);
alpha = CDdatasup(:,2);
CD = CDdatasup(:,3);
% evaluate coefficients
CDcoeffx = [ones(size(CD,1),1), Mn, alpha, alpha.^2] \ CD;
CDexperimental2x = CDcoeffx(1) + CDcoeffx(2).*Mn + CDcoeffx(3).*alpha + ...
    CDcoeffx(4).*alpha.^2; 
CDcoeff2 = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CD;
CDexperimental2 = CDcoeff2(1) + CDcoeff2(2).*Mn + CDcoeff2(3).*alpha + ...
    CDcoeff2(4).*Mn.^2 + CDcoeff2(5).*alpha.^2 + CDcoeff2(6).*Mn.*alpha; 
CDcoeff3 = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
    Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2] \ CD;
CDexperimental3 = CDcoeff3(1) + CDcoeff3(2).*Mn + CDcoeff3(3).*alpha + ...
    CDcoeff3(4).*Mn.^2 + CDcoeff3(5).*alpha.^2 + CDcoeff3(6).*Mn.*alpha + ...
    CDcoeff3(7).*Mn.^3 + CDcoeff3(8).*alpha.^3 + CDcoeff3(9).*Mn.^2.*alpha + CDcoeff3(10).*Mn.*alpha.^2; 
% visualize
figure
hold on
scatter3(Mn, alpha, CD)
% scatter3(Mn, alpha, CDexperimental2)
% scatter3(Mn, alpha, CDexperimental3)
% analyse
R2= corrcoef(CD,CDexperimental2x);
R2(2,1) %.9817
R2= corrcoef(CD,CDexperimental2);
R2(2,1) %.9996
R3 = corrcoef(CD,CDexperimental3);
R3(2,1) %.9999
% visualize result
testpoints = 50;
alpha = linspace(-5,25,testpoints);
Mn = linspace(1,6,testpoints);
CDtest2 = zeros (testpoints,testpoints);
CDtest3=CDtest2;
CDtest2x = CDtest2;
for i = 1:testpoints
    for j = 1:testpoints
        CDtest2x(i,j) = CDcoeffx(1) + CDcoeffx(2).*Mn(j) + CDcoeffx(3).*alpha(i) + ...
            CDcoeffx(4).*alpha(i).^2;
        CDtest2(i,j) = CDcoeff2(1) + CDcoeff2(2).*Mn(j) + CDcoeff2(3).*alpha(i) + ...
            CDcoeff2(4).*Mn(j).^2 + CDcoeff2(5).*alpha(i).^2 +...
            CDcoeff2(6).*Mn(j).*alpha(i);
        CDtest3(i,j) = CDcoeff3(1) + CDcoeff3(2).*Mn(j) + CDcoeff3(3).*alpha(i) + ...
            CDcoeff3(4).*Mn(j).^2 + CDcoeff3(5).*alpha(i).^2 +...
            CDcoeff3(6).*Mn(j).*alpha(i)+ CDcoeff3(7).*Mn(j).^3 + ...
            CDcoeff3(8).*alpha(i).^3 + CDcoeff3(9).*Mn(j).^2.*alpha(i) + CDcoeff3(10).*Mn(j).*alpha(i).^2; 
    end
end
surf(Mn,alpha,CDtest2x,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%}
%{
% CD subsonic
% load xlsx
CDdatasub = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',2 ,'A74:C149');
% purge NaN
CDdatasub(~any(~isnan(CDdatasub), 2),:)=[];
% separate columns
Mn = CDdatasub(:,1);
alpha = CDdatasub(:,2);
CD = CDdatasub(:,3);
% evaluate coefficients
CDcoeffx = [ones(size(CD,1),1), Mn, alpha, alpha.^2] \ CD;
CDexperimental2x = CDcoeffx(1) + CDcoeffx(2).*Mn + CDcoeffx(3).*alpha + ...
    CDcoeffx(4).*alpha.^2; 
CDcoeff2 = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CD;
CDexperimental2 = CDcoeff2(1) + CDcoeff2(2).*Mn + CDcoeff2(3).*alpha + ...
    CDcoeff2(4).*Mn.^2 + CDcoeff2(5).*alpha.^2 + CDcoeff2(6).*Mn.*alpha; 
CDcoeff3 = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
    Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2] \ CD;
CDexperimental3 = CDcoeff3(1) + CDcoeff3(2).*Mn + CDcoeff3(3).*alpha + ...
    CDcoeff3(4).*Mn.^2 + CDcoeff3(5).*alpha.^2 + CDcoeff3(6).*Mn.*alpha + ...
    CDcoeff3(7).*Mn.^3 + CDcoeff3(8).*alpha.^3 + CDcoeff3(9).*Mn.^2.*alpha + CDcoeff3(10).*Mn.*alpha.^2; 
% visualize
figure
hold on
scatter3(Mn, alpha, CD)
% scatter3(Mn, alpha, CDexperimental2)
% scatter3(Mn, alpha, CDexperimental3)
% analyse
R2= corrcoef(CD,CDexperimental2x);
R2(2,1) %.9967
R2= corrcoef(CD,CDexperimental2);
R2(2,1) %.9970
R3 = corrcoef(CD,CDexperimental3);
R3(2,1) %.9985
% visualize result
testpoints = 50;
alpha = linspace(-5,25,testpoints);
Mn = linspace(0,1,testpoints);
CDtest2 = zeros (testpoints,testpoints);
CDtest3=CDtest2;
CDtest2x = CDtest2;
for i = 1:testpoints
    for j = 1:testpoints
        CDtest2x(i,j) = CDcoeffx(1) + CDcoeffx(2).*Mn(j) + CDcoeffx(3).*alpha(i) + ...
            CDcoeffx(4).*alpha(i).^2;
        CDtest2(i,j) = CDcoeff2(1) + CDcoeff2(2).*Mn(j) + CDcoeff2(3).*alpha(i) + ...
            CDcoeff2(4).*Mn(j).^2 + CDcoeff2(5).*alpha(i).^2 +...
            CDcoeff2(6).*Mn(j).*alpha(i);
        CDtest3(i,j) = CDcoeff3(1) + CDcoeff3(2).*Mn(j) + CDcoeff3(3).*alpha(i) + ...
            CDcoeff3(4).*Mn(j).^2 + CDcoeff3(5).*alpha(i).^2 +...
            CDcoeff3(6).*Mn(j).*alpha(i)+ CDcoeff3(7).*Mn(j).^3 + ...
            CDcoeff3(8).*alpha(i).^3 + CDcoeff3(9).*Mn(j).^2.*alpha(i) + CDcoeff3(10).*Mn(j).*alpha(i).^2; 
    end
end
surf(Mn,alpha,CDtest2,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
%}
%{
if Mn > 2
    CDcoeffx = [0.140188997486753,-0.0200905096814203,-0.00200739322354555,0.000621309358144361];
	CD = CDcoeffx(1) + CDcoeffx(2).*Mn + CDcoeffx(3).*alpha + CDcoeffx(4).*alpha.^2;
elseif Mn >1 && Mn<=2
    CDcoeffx = [0.231304857138829,-0.0923845499728899,0.000355575895556675,0.000786190803751548];
	CD = CDcoeffx(1) + CDcoeffx(2).*Mn + CDcoeffx(3).*alpha + CDcoeffx(4).*alpha.^2;    
elseif Mn<=1
    CDcoeffx = [-0.00771662630772907,0.0670075085205795,-0.000381444538875274,0.000919603463455831];
	CD = CDcoeffx(1) + CDcoeffx(2).*Mn + CDcoeffx(3).*alpha + CDcoeffx(4).*alpha.^2;    
end
%}

%% w/o sub-super separation (negative CD)
%CL
%{
% load xlsx
CLdata = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',1 ,'A4:C98');
% purge NaN
CLdata(~any(~isnan(CLdata), 2),:)=[];
% separate columns
Mn = CLdata(:,1);
alpha = CLdata(:,2);
CL = CLdata(:,3);
% evaluate coefficients
CLcoeff = [ones(size(CL,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CL;
CLexperimental2 = CLcoeff(1) + CLcoeff(2).*Mn + CLcoeff(3).*alpha + ...
    CLcoeff(4).*Mn.^2 + CLcoeff(5).*alpha.^2 + CLcoeff(6).*Mn.*alpha; 
CLcoeff = [ones(size(CL,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
    Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2] \ CL;
CLexperimental3 = CLcoeff(1) + CLcoeff(2).*Mn + CLcoeff(3).*alpha + ...
    CLcoeff(4).*Mn.^2 + CLcoeff(5).*alpha.^2 + CLcoeff(6).*Mn.*alpha + ...
    CLcoeff(7).*Mn.^3 + CLcoeff(8).*alpha.^3 + CLcoeff(9).*Mn.^2.*alpha + CLcoeff(10).*Mn.*alpha.^2; 
% visualize
figure
hold on
scatter3(Mn, alpha, CL)
% scatter3(Mn, alpha, CLexperimental2)
% scatter3(Mn, alpha, CLexperimental3)
% analyse
R2= corrcoef(CL,CLexperimental2);
R2(2,1) % 0.9879
R3 = corrcoef(CL,CLexperimental3);
R3(2,1) % 0.9957

testpoints = 50;
CLcoeff = [ones(size(CL,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CL;
alpha = linspace(-5,25,testpoints);
Mn = linspace(0,6,testpoints);
CLtest = zeros (testpoints,testpoints);
for i = 1:testpoints
    for j = 1:testpoints
        CLtest(i,j) = CLcoeff(1) + CLcoeff(2).*Mn(j) + CLcoeff(3).*alpha(i) + ...
            CLcoeff(4).*Mn(j).^2 + CLcoeff(5).*alpha(i).^2 +...
            CLcoeff(6).*Mn(j).*alpha(i);
    end
end
surf(Mn,alpha,CLtest,'EdgeColor','none','LineStyle','none','FaceLighting','phong')


%}
%{
% CLcoeff = [-0.0636596722743266, 0.483834747248333, 0.0643396066871600,...
%     -0.302957181869748,-0.000248775437430644,-0.0136834108337241,...
%     0.0370348358503254,-4.37833607105373e-06,0.000985319452897862,0.000132511094945852];
% 
% CL = CLcoeff(1) + CLcoeff(2).*Mn + CLcoeff(3).*alpha + ...
%     CLcoeff(4).*Mn.^2 + CLcoeff(5).*alpha.^2 + CLcoeff(6).*Mn.*alpha + ...
%     CLcoeff(7).*Mn.^3 + CLcoeff(8).*alpha.^3 + CLcoeff(9).*Mn.^2.*alpha + CLcoeff(10).*Mn.*alpha.^2;
%}
% CD
%{
% load xlsx
CDdata = xlsread('X34_aero (ref Brauckmann 1999 JSR).xlsx',2 ,'A3:C149');
% purge NaN
CDdata(~any(~isnan(CDdata), 2),:)=[];
% separate columns
Mn = CDdata(:,1);
alpha = CDdata(:,2);
CD = CDdata(:,3);
% evaluate coefficients
CDcoeff = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CD;
CDexperimental2 = CDcoeff(1) + CDcoeff(2).*Mn + CDcoeff(3).*alpha + ...
    CDcoeff(4).*Mn.^2 + CDcoeff(5).*alpha.^2 + CDcoeff(6).*Mn.*alpha; 
CDcoeff = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
    Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2] \ CD;
CDexperimental3 = CDcoeff(1) + CDcoeff(2).*Mn + CDcoeff(3).*alpha + ...
    CDcoeff(4).*Mn.^2 + CDcoeff(5).*alpha.^2 + CDcoeff(6).*Mn.*alpha + ...
    CDcoeff(7).*Mn.^3 + CDcoeff(8).*alpha.^3 + CDcoeff(9).*Mn.^2.*alpha + CDcoeff(10).*Mn.*alpha.^2; 
CDcoeff = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha,...
    Mn.^3, alpha.^3, Mn.^2.*alpha, Mn.*alpha.^2,...
    Mn.^4, alpha.^4, Mn.^3.*alpha, Mn.*alpha.^3, Mn.^2.*alpha.^2] \ CD;
CDexperimental4 =CDcoeff(1) + CDcoeff(2).*Mn + CDcoeff(3).*alpha + ...
    CDcoeff(4).*Mn.^2 + CDcoeff(5).*alpha.^2 + CDcoeff(6).*Mn.*alpha + ...
    CDcoeff(7).*Mn.^3 + CDcoeff(8).*alpha.^3 + CDcoeff(9).*Mn.^2.*alpha + CDcoeff(10).*Mn.*alpha.^2 + ...
    CDcoeff(11).*Mn.^4 + CDcoeff(12).*alpha.^4 + CDcoeff(13).*Mn.^3.*alpha + ...
    CDcoeff(14).*Mn.*alpha.^3 +  CDcoeff(15).*Mn.^2.*alpha.^2; 
% visualize
figure
hold on
scatter3(Mn, alpha, CD)
% scatter3(Mn, alpha, CDexperimental2)
% scatter3(Mn, alpha, CDexperimental3)
% scatter3(Mn, alpha, CDexperimental4)
% analyse
R2= corrcoef(CD,CDexperimental2);
R2(2,1) % 0.9770
R3 = corrcoef(CD,CDexperimental3);
R3(2,1) % 0.9881
R4 = corrcoef(CD,CDexperimental4);
R4(2,1) % 0.9980


testpoints = 50;
CDcoeff = [ones(size(CD,1),1), Mn, alpha, Mn.^2, alpha.^2, Mn.*alpha] \ CD;
alpha = linspace(-5,25,testpoints);
Mn = linspace(0,6,testpoints);
CDtest = zeros (testpoints,testpoints);
for i = 1:testpoints
    for j = 1:testpoints
        CDtest(i,j) = CDcoeff(1) + CDcoeff(2).*Mn(j) + CDcoeff(3).*alpha(i) + ...
            CDcoeff(4).*Mn(j).^2 + CDcoeff(5).*alpha(i).^2 +...
            CDcoeff(6).*Mn(j).*alpha(i);
    end
end
surf(Mn,alpha,CDtest,'EdgeColor','none','LineStyle','none','FaceLighting','phong')

%}
%{
% CDcoeff = [0.413058140353296,-1.68853504423995,-0.00877072497428170,2.16624753778642,...
% 0.00104357462389481,0.0190918630988641,-0.896360294715539,9.73780385415617e-06,...
% -0.00922482218215578,-0.000391375623782086,0.0967599325572793,-2.54805046366662e-07,...
% 0.00103446330194631,1.76615243708157e-06,4.26334462345704e-05];
% 
% CD = CDcoeff(1) + CDcoeff(2).*Mn + CDcoeff(3).*alpha + ...
%     CDcoeff(4).*Mn.^2 + CDcoeff(5).*alpha.^2 + CDcoeff(6).*Mn.*alpha + ...
%     CDcoeff(7).*Mn.^3 + CDcoeff(8).*alpha.^3 + CDcoeff(9).*Mn.^2.*alpha + CDcoeff(10).*Mn.*alpha.^2 + ...
%     CDcoeff(11).*Mn.^4 + CDcoeff(12).*alpha.^4 + CDcoeff(13).*Mn.^3.*alpha + ...
%     CDcoeff(14).*Mn.*alpha.^3 +  CDcoeff(15).*Mn.^2.*alpha.^2;
%}