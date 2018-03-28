function [res] = extract_control_law (param, c_opt, phases)
%{
Function to extract the control law from the optimisation vector
INPUT:
-finaltime: xf(4), target flight time
-param: parameters characterizing the problem
-c_opt: the optimisation vector

OUTPUT:
res(i).t: each phase's linspaced time vector for the integration
res(i).control: each phase's control matrix, first column is linspaced
time, second to last are controls (eg: alpha, throttle and bank are 3 columns)

(c) 2015, F Toso, Centre for Future Air-Space Transportation
Technology, Univeristy of Strathclyde
%}
k=1;
index = 1;
for ip = 1 : param.np                                                       % select phase
    phases(ip).odet = phases(ip).nc*2;                                        % odet moved from  
    for ine = 1 : phases(ip).ne                                             % select element
        %% extract res.t and res.control
        res(index).control=zeros(phases(ip).nc, param.nv+1);                % initialize control matrix 
%         extract control values from opt vector
        for iv = 1 : phases(ip).nu
            if phases(ip).cbounds(iv,1)~=phases(ip).cbounds(iv,2)
%                 fill the row with the optimised values, advance k
                res(index).control(:,1+iv) = c_opt(k:k-1+phases(ip).nc);
                k=k+phases(ip).nc;
            else
%                 if not optimizable, fill the column with the first guess
                res(index).control(:,1+iv) = phases(ip).cbounds(iv,1)*ones(phases(ip).nc,1);
            end
        end
%         extract time value from opt vector
        tof = max(0.1,c_opt(k));                                                     % calc time of flight, fix for out of bound checks (t<0)
        tsegments = max([3,phases(ip).odet,ceil(tof/phases(ip).tstep)]);        % evaluate the number of points for integration
        res(index).t=linspace(0,tof,tsegments);                             % linspaced time vector for integration
        res(index).control(:,1)=phases(ip).dist(0, tof, phases(ip).nc)';    % create first column of the control matrix with time
        k=k+1;
        %% extract res.x0
        if isempty(phases(ip).x0_fixed) || ine ~=1                     % if all x0 are optimisation variables
            if ip>1 && ine ==1                                          % staging
                if phases(ip).vehicle.gtow==phases(ip-1).vehicle.gtow   % no staging, optimise mass
                    res(index).x0=c_opt(k:k+6);                         % matching
                    k=k+7;                                              % matching
                else                                                    % if there is staging, mass is constrained
                    res(index).x0=[c_opt(k:k+5), phases(ip).vehicle.gtow];% staging
                    k=k+6;                                              % staging
                    if isfield(param,'Scale_S2_GTOW') && ~exist('gtow_s2_opt','var') 
                        if param.Scale_S2_GTOW(1)-param.Scale_S2_GTOW(2)~=0
                            res(index).x0(7) = c_opt(k);
                            k=k+1;
                            gtow_s2_opt = 1;
                        else
                            res(index).x0(7)=param.Scale_S2_GTOW(1);
                        end
                    end
                end
            else
                res(index).x0=c_opt(k:k+6);                             % matching
                k=k+7;                                                  % matching
            end
        else
            for ins = 1 : size(phases(ip).xbounds,1)                    % for each state variable
                if phases(ip).x0_fixed(ins)~=1                         % if it is not fixed
                        res(index).x0(ins)=c_opt(k);                    % if it's not fixed, take value from c_opt
                        k=k+1;
                else
                    res(index).x0(ins)=phases(ip).x0(ins);
                end
            end
        end
        index = index+1;
    end
end
