function [extracted] = interp_cos (control,time)

%{
Function to interpolate the control matrix

INPUT
control = nc*(1+nv); control law vector, first column is time
time = sample time 

(c) 2015, F Toso, Centre for Future Air-Space Transportation Technology, Univeristy of Strathclyde
%}

% cosine interpolation

tindex=1;
while time>=control(tindex,1) && tindex < size(control,1)
    tindex=tindex+1;
end

percentage=(time-control(tindex-1,1))/(control(tindex,1)-control(tindex-1,1));
cosper = (1-cos(percentage*pi))/2;

extracted = zeros(1,size(control,2)-1);

for i = 2:size(control,2)
    extracted(i-1) = control(tindex-1,i)+cosper*(control(tindex,i)-control(tindex-1,i));
end