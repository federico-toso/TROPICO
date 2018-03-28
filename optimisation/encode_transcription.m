function [LB, UB, phases,param] = generate_bounds (param, phases)
%{
function to create the optimisation vector.

INPUT:
- c0 }
- LB }} Should be == []
- UB }
- param
- phases
OUTPUT:
- c0 created
- LB created
- UB created
- phases

(c) 2016, F Toso, R Garner, Centre for Future Air-Space Transportation
Technology, Univeristy of Strathclyde
%}

%% Initialisation

LB=[];
UB=[];
c0=[];
param.np = numel(phases);

%% Transcription of state, control & time

for ip = 1 : param.np                                                       % select phase
    phases(ip).nu = size(phases(ip).cbounds,1);                             % Number of elements in the control vector, e.g., [angle of attack, throttle, bank]
    for ine = 1 : phases(ip).ne                                             % select element
%         add control values, and bounds
        for iv = 1:phases(ip).nu
            if phases(ip).cbounds(iv,1)~=phases(ip).cbounds(iv,2)
%                 c0 = [c0, phases(ip).fg(iv)         *ones(1,phases(ip).nc)];    % first guess
                LB = [LB, phases(ip).cbounds(iv,1)  *ones(1,phases(ip).nc)];    % lower bound
                UB = [UB, phases(ip).cbounds(iv,2)  *ones(1,phases(ip).nc)];    % upper bound
            end
        end
%         add time variable
        phases(ip).tof(1) = max(phases(ip).tof(1),1);
%         c0 = [c0, (phases(ip).tof(1)+phases(ip).tof(2))/2/phases(ip).ne];   %first guess is average between bounds
        LB = [LB, phases(ip).tof(1)/phases(ip).ne];
        UB = [UB, phases(ip).tof(2)/phases(ip).ne];
%         add state variables
        if isempty(phases(ip).x0_no_opt) || ine ~=1                         % if x0_no_opt is not defined, add all
            LB = [LB, phases(ip).xbounds(:, 1)'];
            UB = [UB, phases(ip).xbounds(:, 2)'];
%             c0 = [c0, (phases(ip).xbounds(:, 2)+phases(ip).xbounds(:, 1))'/2];
        else                                                                % if not, add the specified ones
            for ins = 1 : size(phases(ip).xbounds,1)                        % for each state variable
                if phases(ip).x0_no_opt(ins)~=1 
%                     c0 = [c0, phases(ip).x0(ins)];                          % if it's not fixed, add value
                    LB = [LB, phases(ip).xbounds(ins, 1)];
                    UB = [UB, phases(ip).xbounds(ins, 2)];
                end
            end
        end
        if ip>1 && ine ==1                                                  % in case of staging, don't optimise mass, remove last element
            if phases(ip).vehicle.gtow~=phases(ip-1).vehicle.gtow
%                 c0=c0(1:end-1);   
                LB=LB(1:end-1);    UB=UB(1:end-1);
            	if isfield(param,'Scale_S2_GTOW') && ~exist('gtow_s2_opt','var')   % second stage is optimized, add values for it
                    if param.Scale_S2_GTOW(1)-param.Scale_S2_GTOW(2)~=0             % if optimisation bounds are different
                        gtow_s2_opt=1;
%                         c0 = [c0, mean(param.Scale_S2_GTOW)];
                        LB = [LB, param.Scale_S2_GTOW(1)];
                        UB = [UB, param.Scale_S2_GTOW(2)];
%                         param.Position_S2_GTOW = length(UB);
                        if param.Scale_S2_GTOW(2)>phases(ip).vehicle.gtow || param.Scale_S2_GTOW(1)<phases(ip).vehicle.m0
                            error(['phases(',num2str(ip),...
                                ') xbounds and vehicle bounds are more strict then S2_GTOW ones. Usually the solution lies in changing gtow value in "load_vehicle.m" to accomodate for the "param.Scale_S2_GTOW" bounds'])
                        end
                    else
                        disp('WARNING: S2_GTOW bounds have the same value and will overwrite vehicle.gtow. To avoid, comment out the definition of "param.Scale_S2_GTOW"')
                    end
                end
            end
        end
    end
end

%% Transcription of additional extra optimisation variables for scaling of forces (L,D, T, mp)

if isfield(param,'Scale_S1_Engine') || isfield(param,'Scale_S2_Engine') || isfield(param,'Scale_S1_Surf') ||...
        isfield(param,'Scale_S2_Surf') || isfield(param,'Scale_S2_GTOW')
    param.ScalingActive=1;
    [LB, UB, phases] = AddScalingPrameters (LB, UB, param, phases);
else
    param.ScalingActive=0;
end

%% Scale the opt vector between 0 & 1            

param.LB=LB;
param.UB=UB;
% c0=(c0-LB)./(UB-LB);
UB=ones(1,numel(UB));
LB=0*UB;


