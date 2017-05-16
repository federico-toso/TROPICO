function [extracted] = interp_pch (control,time)
%#codegen
%{
[extracted] = interp_pch (control,time)

Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) interpolation of a
control matrix. 

INPUT
control = nx4 control matrix, 1st column is the time grid, next 3 columns
are values at each time point
time = point where the interpolation has to be done

OUTPUT
extracted = vector containing the interpolated points

(c) 2015, F Toso, Centre for Future Air-Space Transportation Technology, Univeristy of Strathclyde
%}

% if first point
if time == 0
    extracted = control(1,2:end);
    return
end

tindex=1;
% find next datapoint
while time>=control(tindex,1) && tindex < size(control,1)
    tindex=tindex+1;
end

step = control(tindex,1)-control(tindex-1,1);

% initialize the interpolated values vector and fill it with pchip method
extracted = zeros(1,size(control,2)-1);
% extractedcheck = extracted;
for i = 2:size(control,2)

    x1=control(tindex-1,1);
    x2=control(tindex,1);
    y1=control(tindex-1,i);
    y2=control(tindex,i);

    d12 = (y2-y1)/step;

    if tindex == 2
        d01 = 0;
        d23 = (control(tindex+1,i) - y2)/step;
    elseif tindex == size(control,1)
        d01 = (y1-control(tindex-2,i))/step;
        d23 = 0;
    else
        d01 = (y1-control(tindex-2,i))/step;
        d23 = (control(tindex+1,i) - y2)/step;
    end
    
    if sign(d01)==sign(d12)     % monotonic function
        d1 = (d01+d12)/2;
    else
        d1 = 0;
    end
% 
    if sign(d12)==sign(d23)     % monotonic function
        d2 = (d12+d23)/2;
    else
        d2 = 0;
    end

% this method is faster than matlab's, valid for a 4 column control matrix
    den = (x1 - x2)*(x1^2 - 2*x1*x2 + x2^2);
    a = -(2*y1 - 2*y2 - d1*x1 + d1*x2 - d2*x1 + d2*x2)/den;
    b = (3*x1*y1 - 3*x1*y2 + 3*x2*y1 - 3*x2*y2 - d1*x1^2 + 2*d1*x2^2 - 2*d2*x1^2 + d2*x2^2 - d1*x1*x2 + d2*x1*x2)/den;
    c = -(d1*x2^3 - d2*x1^3 + d1*x1*x2^2 - 2*d1*x1^2*x2 + 2*d2*x1*x2^2 - d2*x1^2*x2 + 6*x1*x2*y1 - 6*x1*x2*y2)/den;
    d = (x1^3*y2 - x2^3*y1 + d1*x1*x2^3 - d2*x1^3*x2 + 3*x1*x2^2*y1 - 3*x1^2*x2*y2 - d1*x1^2*x2^2 + d2*x1^2*x2^2)/den;
    
    extracted(i-1) = a*time^3 + b*time^2 + c*time + d;
    
%     A =[x1^3   x1^2 x1 1;
%         3*x1^2 2*x1 1  0;
%         x2^3   x2^2 x2 1;
%         3*x2^2 2*x2 1  0];
%     B =[y1; d1; y2; d2];
%     X = A\B;
%     extracted(i-1) = X(1)*time^3 + X(2)*time^2 + X(3)*time + X(4);
end