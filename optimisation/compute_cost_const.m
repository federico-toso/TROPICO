function [cost_function,c,ceq] = compute_cost_const(u_opt, param, phases)
%{
function to calculate the cost function and the constraint vectors

INPUT:
-c_opt: the optimisation vector to evaluate
-x0: the starting vector, scaled
-xf: target state vector, not scaled
-totime: target flight time
-x_scale: scaling vector
-param: parameters characterizing the problem

OUTPUT:
-cost_function: for the optimisation algorithm
-c: vector of inequality constraints
-ceq: vector of the equality constraints

(c) 2016, F Toso, R Garner, Centre for Future Air-Space Transportation
Technology, Univeristy of Strathclyde
%}

loadconstants
c=[];       ceq=[];         

%% extract control law, integration points(t) and x0
c_opt=param.LB+u_opt.*(param.UB-param.LB);

% optional scaling routine for engine and surfaces
[phases] = ExtractScalingPrameters (param, c_opt, phases);

[res] = extract_control_law (param, c_opt, phases);                          % initialize control matrix, time vector and x0 vector
[res] = propagate_trajectory (param, phases, const, res);

%% param.constraint_debug
if isfield(param,'constraint_debug')
    verbosity = 1;
else
    verbosity = 0;
end

%% control matching between trajectory elements (starts from second element)
index=0;
for ip = 1 : param.np                                                       % select phase
    if ~isempty(phases(ip).control_tol)
        controls_scale = phases(ip).control_tol*1e4;                        % control scaling based on manual input
    else
        controls_scale = phases(ip).cbounds(:,2)-phases(ip).cbounds(:,1);   % default control scaling on the range between bounds
    end
    for ine = 1 : phases(ip).ne                                             % select element
        index = index+1;                                                    
        if index > 1                                                        % after first element, match controls
            if ine == 1 && ~isempty(phases(ip).continue_from)
                link_to=0;                                                  % select the previous element
                for ilink = 1:phases(ip).continue_from
                    link_to = link_to + phases(ilink).ne;
                end
                after_phase = phases(ip).continue_from;
            else
                link_to = index-1;
                after_phase = (ine==1)*(ip-1);
            end
            if ine == 1 && ~isempty(phases(ip).controls_no_match)           % match some of the controls, or
                for ic = 1 : size (phases(ip).cbounds,1)
                    if phases(ip).controls_no_match(ic)~=1
                        if phases(ip).cbounds(ic,2)~=phases(ip).cbounds(ic,1)
                            ceq = [ceq; (res(index).control(1,ic+1)-res(link_to).control(end,ic+1))/controls_scale(ic)];
                        end
                    end
                end
            else                                                            % match all controls
                for ic = 1 : size (phases(ip).cbounds,1)
                    if phases(ip).cbounds(ic,2)~=phases(ip).cbounds(ic,1)
                        ceq = [ceq; (res(index).control(1,ic+1)-res(link_to).control(end,ic+1))/controls_scale(ic)];
                    elseif after_phase~=0
                        if phases(after_phase).cbounds(ic,2)~=phases(after_phase).cbounds(ic,1)
                            ceq = [ceq; (res(index).control(1,ic+1)-res(link_to).control(end,ic+1))/...
                                (phases(after_phase).cbounds(ic,2)-phases(after_phase).cbounds(ic,1))];
                        end
                    end
                end
            end
        end
    end
end

if verbosity == 1
    counter_c(1)=size(c,1);
    counter_ceq(1)=size(ceq,1);
end

%% state matching between trajectory elements (starts from second element, includes staging)
index = 0;
states_scale = zeros(param.ns,param.np);
for ip = 1 : param.np 
    if ~isempty(phases(ip).state_tol)
        states_scale(:,ip) = phases(ip).state_tol/1e4;
    else
        states_scale(:,ip) = max([phases(ip).xbounds(:,2)-phases(ip).xbounds(:,1)],1);         % state scaling
    end
    for ine = 1 : phases(ip).ne
        index = index+1;
        if index > 1 % match the final states of an element with the start of the next one
            if ine == 1 && ~isempty(phases(ip).continue_from)
                link_to=0;
                for ilink = 1:phases(ip).continue_from
                    link_to = link_to + phases(ilink).ne;
                end
            else
                link_to = index-1;
            end
            ceq = [ceq; (res(index).x0'-res(link_to).x(end,:)')./states_scale(:,ip)];
        end
        if ip >1 && ine==1 % remove the mass equality to satisfy staging
            if phases(ip).vehicle.gtow~=phases(ip-1).vehicle.gtow
                ceq = ceq(1:end-1);     % remove constraint on mass matching if staging
                if link_to+1~=index     %add staging constraint
                    k=length(ceq);
                    k=k+1-7*(index-link_to-1);
                    ceq(k)=(res(index).x0(7)+res(link_to+1).x0(7))/res(link_to).x(end,7)-1;
                end
            end
        end
    end
end

if verbosity == 1
    counter_c(2)=size(c,1);
    counter_ceq(2)=size(ceq,1);
end
%% constrain vehicle final mass > m0 (dry mass)
index = 0;
for ip = 1:param.np
    index = index+phases(ip).ne;
    for ip2 = ip:param.np
        if isequal(phases(ip).vehicle,phases(ip2).vehicle) && ip~=ip2
            break
        end
        if ip2 == param.np
            c = [c; (phases(ip).vehicle.m0-res(index).x(end,7))/phases(ip).vehicle.gtow];
        end
    end
end

if verbosity == 1
    counter_c(3)=size(c,1);
    counter_ceq(3)=size(ceq,1);
end
%% end phase matching (xf values)
index=0;
for ip = 1 : param.np                                                       % select phase
    index=index+phases(ip).ne;
    if ~isempty(phases(ip).xf_fixed)
        for ifinal = 1 : length(phases(ip).xf_fixed)
            switch phases(ip).xf_fixed(ifinal)
            case 1  % equality
                ceq = [ceq; (res(index).x(end,ifinal)-phases(ip).xf(ifinal))/states_scale(ifinal,ip)];
            case 2  % more than
                c = [c; (phases(ip).xf(ifinal)-res(index).x(end,ifinal))/states_scale(ifinal,ip)];
            case -1 % less than
                c = [c; (-phases(ip).xf(ifinal)+res(index).x(end,ifinal))/states_scale(ifinal,ip)];
            end
        end
    end
    if ~isempty(phases(ip).xf_range)
        for ibound = 1:param.ns
            if phases(ip).xf_range(ibound,1)~=phases(ip).xf_range(ibound,2)
                lower_diff=(phases(ip).xf_range(ibound,1)-res(index).x(end,ibound));
                upper_diff=(res(index).x(end,ibound)-phases(ip).xf_range(ibound,2));
                c = [c; max([0, lower_diff, upper_diff])/states_scale(ibound,ip)];
            end
        end
    end
 
end

if verbosity == 1
    counter_c(4)=size(c,1);
    counter_ceq(4)=size(ceq,1);
end
%% evaluate states, constraints, ecc.
accelerations=[];    thermal=[];  forces =[];   controls =[]; Mn =[]; atmo =[];
index=0;
for ip = 1 : param.np                                                       % select phase
    for ine = 1 : phases(ip).ne                                             % select element
        index = index+1;                                                    % propagate trajectory
        for iicc=1:size(res(index).x,1)      
            [~, con]= phases(ip).dynamics(res(index).t(iicc), res(index).x(iicc,:), res(index).control, const, phases(ip));
            accelerations=[accelerations; con.acc]; %[accx, accz]
            thermal=[thermal; con.temp];            %[temperature, heat]
            forces=[forces; con.forces];            %[L,D,FT,q]
            controls=[controls; con.controls];
            Mn=[Mn; con.Mn];
            atmo = [atmo; con.atmo];
        end
    end
end

%% constraints on the whole trajectory
if isfield(param,'c')
    for ic = 1:size(param.c,1)
        switch param.c{ic,1}
        case 'c_acc' % constrain accelerations
            scaleacc = 1e0*const.g0;
            switch length(param.c{ic,2})
            case 1 % if only one value is set, vector sum is computed
                c = [c; (max(hypot(accelerations(:,1),accelerations(:,2)))-param.c{ic,2}*const.g0)/scaleacc]; % i.e. {'c_acc',6}
            case 2 % if two values, first is alon accx, second is accz
                c = [c; (max(abs(accelerations(:,1))-param.c{ic,2}(1)*const.g0))/scaleacc; % i.e. {'c_acc',[6 6]}
                        (max(abs(accelerations(:,2))-param.c{ic,2}(2)*const.g0))/scaleacc];
            end
        case 'c_temp' % constrain temperature
            c = [c; max(thermal)/param.c{ic,2}(1)-1];                           % i.e. {'c_temp',3000}
        case 'c_maxq' % constrain max dynamic pressure
            c = [c; max(forces(:,4))/param.c{ic,2}(1)-1];                         % i.e. {'c_maxq',2e6}
        end
    end
end

if verbosity == 1
    counter_c(5)=size(c,1);
    counter_ceq(5)=size(ceq,1);
end
%% constraints on the single phases
index = 0;
for ip = 1:param.np
    index = index + phases(ip).ne;
    if isfield (phases(ip),'ceq')
        for ic = 1:size(phases(ip).ceq,1)
            [a,e,i] = enu2orbital (res(index).x(end,:)); %res(index).a,res(index).e,res(index).i
            switch cell2mat(phases(ip).ceq(ic,1))
            case 'ceq_circ' % constrain eccentricity
                e = min(max(0,e),e^(1/5));
                ceq = [ceq; e-phases(ip).ceq{ic,2}];        % i.e. {'ceq_e',0}
            case 'ceq_i' % constrain incliantion
                ceq = [ceq; (i-phases(ip).ceq{ic,2})/pi];     % i.e. {'ceq_i',deg2rad(88)}
            case 'ceq_sso' % constrain sun synchronous
                a = max(0,a);
                ceq = [ceq; cos(i)+(a/12352e3)^(7/2)];      % i.e. {'ceq_sso',[]}
            end
        end
    end
end

if verbosity == 1
    counter_c(6)=size(c,1);
    counter_ceq(6)=size(ceq,1);
end

%% Cost Function Calculation
switch param.objfun
case 'max_mass_frac'     % maximize final vehicle mass FRACTION
    cost_function = -res(end).x(end,7)/res(1).x(1,7);
case 'max_payload'   % maximize final mass
    cost_function = -res(end).x(end,7)/1e3;
case 'min_gtow'      % minimize initial mass
    cost_function = res(1).x(1,7)/1e5;
case 'max_altitude'     % maximize final altitude
	cost_function = -res(end).x(end,1)/1e6; 
case 'max_speed'     % maximize final speed
	cost_function = -res(end).x(end,2)/1e4;   
case 'min_peak_heat' % minimize max temperature
    cost_function = max(thermal)/3000;
case 'min_max_acc'   % minimize max accelerations
    cost_function = max([max(accelerations(:,1)),max(accelerations(:,2))])/const.g0;
case 'max_c3'        % maximize final vehicle characteristic energy
    cost_function = (0.5*res(end).x(end,2)^2-const.mu/(res(end).x(end,1)+const.rE))/...
                 (0.5*res(1).x(1,2)^2-const.mu/(res(1).x(1,1)+const.rE));
case 'min_fuel'      % minimise fuel/stage fraction
    cost_function = 0;
    for icost = 1 : numel(res)
        cost_function = cost_function-(res(icost).x(end,7)-res(icost).x(1,7))/res(icost).x(1,7);
    end
case 'min_time'      % minimize total time
    cost_function = 0;
    for icost=1:numel(res), cost_function = cost_function + res(icost).t(end)/60; end
case 'max_distance'      % maximize distance
    cost_function = -distance_haversine(res(1).x(1,5), res(1).x(1,6),...
        res(end).x(end,5), res(end).x(end,6), const.rE)/const.rE;
case 'min_distance'      % minimize distance
    cost_function = distance_haversine(res(1).x(1,5), res(1).x(1,6),...
        res(end).x(end,5), res(end).x(end,6), const.rE)/const.rE;
case 'max_latitude'   % maximize final latitude
    cost_function = -res(end).x(end,5);
case 'feasibility'   % minimize constraints
    cost_function = sum(c(c>0)) + sum(abs(ceq(ceq~=0)));
    c=[]; ceq=[];
case 'custom'   % custom cost max range min weight
%     cost_function = -res(1).x(end,1)/1e5; % max separation altitude
%     cost_function1 = res(1).x(1,7)/1e5; % min gtow
    cost_function = -distance_haversine(res(1).x(1,5), res(1).x(1,6),...
        res(4).x(end,5), res(4).x(end,6), const.rE)/const.rE; % max distance on res(end-3)
%     cost_function = param.scalecost*cost_function1 + (1-param.scalecost)*cost_function2;
otherwise
    error('No cost function defined, review input settings')
end

%% path constraint, no fly zones
% param.no_fly contains lat lon of the centre of the circle, radius
if isfield(param,'no_fly')
    if ~isempty(param.no_fly)
        scale_distance = 1e6;
        distance=[];
        for index_nfz=1:numel(res)                  % first, evaluate distances
            segment_distance = distance_haversine(res(index_nfz).x(:,5), res(index_nfz).x(:,6),...
                                                     param.no_fly(:,1), param.no_fly(:,2), const.rE);
            for nfz = 1:size(param.no_fly,1)        % then, replace the values above altitude treshold with safety margins
                segment_distance(res(index_nfz).x(:,1)>param.no_fly(nfz,4),nfz) = param.no_fly(nfz,3)+1e5; % [possible update includes doing it only if inside the circle]
            end
            distance=[distance;segment_distance];   % include eache segment in a general matrix
        end
        min_distance = min(distance)/scale_distance;% calculate minimum distances
        c=[c;param.no_fly(:,3)/scale_distance-min_distance'];
    end
end

%% final checks for nans
if sum(isnan(c(:))) >0 || sum(isnan(ceq(:))) >0 || sum(isnan(cost_function)) >0
    disp('nan in compute_cost_const, check end of file')
    keyboard
elseif sum(isinf(c(:))) >0 || sum(isinf(ceq(:))) >0 || sum(isinf(cost_function)) >0
    disp ('inf in compute_cost_const, check end of file')
    keyboard
end

%% debug info
if verbosity == 1
    debugcell_c = cell(size(c,1),3);
    debugcell_ceq = cell(size(c,1),3);
    for loopi =  1:size(c,1)
        debugcell_c{loopi,1} = loopi;
        debugcell_c{loopi,3} = c(loopi);
        if sum(loopi==counter_c)>0
            section_idx = find(counter_c==loopi,1,'first');
            if section_idx==1
                debug_string='controls';
            elseif section_idx==2
                debug_string='states';
            elseif section_idx==3
                debug_string='positive mass';
            elseif section_idx==4
                debug_string='xf';
            elseif section_idx==5
                debug_string='whole traj constraints';
            elseif section_idx==6
                debug_string='single phase constraints';
            end
            debugcell_c{loopi,2} = debug_string;
        end
%         if debugcell_c{loopi,3}>1e-4
%             debugcell_c{loopi,4} = debugcell_c{loopi,3};
%         else
%             debugcell_c{loopi,4}=[];
%         end
    end
    for loopi =  1:size(ceq,1)
        debugcell_ceq{loopi,1} = loopi;
        debugcell_ceq{loopi,3} = ceq(loopi);
        if sum(loopi==counter_ceq)>0
            section_idx = find(counter_ceq==loopi,1,'first');
            if section_idx==1
                debug_string='controls';
            elseif section_idx==2
                debug_string='states';
            elseif section_idx==3
                debug_string='positive mass';
            elseif section_idx==4
                debug_string='xf';
            elseif section_idx==5
                debug_string='whole traj constraints';
            elseif section_idx==6
                debug_string='single phase constraints';
            end
            debugcell_ceq{loopi,2} = debug_string;
        end
%         if abs(debugcell_ceq{loopi,3})>1e-4
%             debugcell_ceq{loopi,4} = abs(debugcell_ceq{loopi,3});
%         else
%             debugcell_ceq{loopi,4}=[];
%         end
    end
    disp('c:')
    disp(debugcell_c)
    disp('ceq:')
    disp(debugcell_ceq)
    debugcell_cost = {param.objfun,cost_function};
    disp('cost:')
    disp(debugcell_cost);
end

%% Debug
%  c,ceq, cost_function % keyboard
% if isempty(ceq) || isempty(c)
%     keyboard
% end
% if max(isinf(ceq)) || max(isinf(c))
%     keyboard
% end


%% AIAA HYP 2017 specific constraints
% % ceq = [ceq; cos(i)]; % polar
% airport_distance = distance_haversine(res(3).x(end,5), res(3).x(end,6),...
%                    param.landing_zones(:,1), param.landing_zones(:,2), const.rE);
% min_airport_dist=(min(airport_distance));
% % next fuinction has value 0 near the airport
% sigmoid_distance = min_airport_dist/100000 * 1/(1+(e^(-1/100*(min_airport_dist-10000))));
% ceq = [ceq; sigmoid_distance];

%% IAC2016 specific constraints
% m0_325t=0.15;
% m0_100t=0.25;
% slope = (m0_325t-m0_100t)/(325-100);
% coeff = m0_325t - slope*325;
% min_mass_fraction = slope * res(1).x(1,7)/1e3 + coeff;
% index=0;                                                % find last segment of S1 flight
% for ip = 1 : param.np
%     for ine = 1 : phases(ip).ne
%         index = index+(phases(ip).vehicle.Sgross==405);
%     end
% end
% % % % % % % c = [c; (min_mass_factor*res(1).x(1,7)-res(end).x(end,7))/res(1).x(1,7)]; %final mass must be 20% or more of the initial one
% % c = [c; min_mass_fraction-(res(index).x(end,7)-1e3)/res(1).x(1,7)]; %final mass must be 2X% or more of the initial one
% c = [c; min_mass_fraction-(res(4).x(end,7)-1e3)/res(1).x(1,7)]; %final mass must be 2X% or more of the initial one
