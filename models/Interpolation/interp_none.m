function [extracted] = interp_none(control,~)

%{
Function to interpolate the control matrix

INPUT
control = nc*(1+nv); control law vector, first column is time
time = sample time 

(c) 2015, F Toso, Centre for Future Air-Space Transportation Technology, Univeristy of Strathclyde
%}

% no interpolation, the first row of the control matrix is extracted
extracted=control(1,2:end);