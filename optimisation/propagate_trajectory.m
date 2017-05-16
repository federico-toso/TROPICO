function [res] = propagate_trajectory (param, phases, const, res)

index=0;
for ip = 1 : param.np                                                       % select phase
    for ine = 1 : phases(ip).ne                                             % select element
        index = index+1;                                                    % propagate trajectory
        if isequal(phases(ip).ode,@ode1) || isequal(phases(ip).ode,@ode2) || isequal(phases(ip).ode,@ode3) ||isequal(phases(ip).ode,@ode4) ||isequal(phases(ip).ode,@ode5)
            res(index).x=phases(ip).ode(phases(ip).dynamics,res(index).t, res(index).x0, res(index).control, const, phases(ip));
        else %use matlab's ode
            [~,res(index).x]=phases(ip).ode(phases(ip).dynamics,res(index).t, res(index).x0, [], res(index).control, const, phases(ip));
        end
        if size(res(index).t,2) ~= size(res(index).x,1) %this routine fixes the output of a stiff problem solved with matlab's ode
            res(index).x=[res(index).x;...
                repmat(res(index).x(end,:),size(res(index).t,2)-size(res(index).x,1),1)];
        end
        if index > 1 %&& res(index).t(1)==0
            if ine == 1 && ~isempty(phases(ip).continue_from)
                link_to=0;
                for ilink = 1:phases(ip).continue_from
                    link_to = link_to + phases(ilink).ne;
                end
            else
                link_to = index-1;
            end
            res(index).t = res(index).t+res(link_to).t(end);
            res(index).control(:,1) = res(index).control(:,1)+res(link_to).t(end);
        end
    end
end

