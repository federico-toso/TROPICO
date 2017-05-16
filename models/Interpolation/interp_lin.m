function [extracted] = interp_lin(control,time)

%{
Function to interpolate the control matrix

INPUT
control = nc*(1+nv); control law vector, first column is time
time = sample time 

(c) 2015, F Toso, Centre for Future Air-Space Transportation Technology, Univeristy of Strathclyde
%}

% linear interpolation

tindex=1;
while time>=control(tindex,1) && tindex < size(control,1)
    tindex=tindex+1;
end

extracted = zeros(1,size(control,2)-1);
percentage=(time-control(tindex-1,1))/(control(tindex,1)-control(tindex-1,1));

for i = 2:size(control,2)
    extracted(i-1)=control(tindex-1,i)+percentage*(control(tindex,i)-control(tindex-1,i));
%     alpha=control(tindex-1,2)+percentage*(control(tindex,2)-control(tindex-1,2));
%     throttle=control(tindex-1,3)+percentage*(control(tindex,3)-control(tindex-1,3));
%     bank=control(tindex-1,4)+percentage*(control(tindex,4)-control(tindex-1,4));
end