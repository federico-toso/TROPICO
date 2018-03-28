function eval_plot (c0, param, phases)
%{
function to plot a optimisation vector

INPUT:
-c0: the optimisation vector to evaluate
-x0: the starting point, not scaled
-totime: target flight time
-x_scale: scaling vector
-param: parameters characterizing the problem

OUTPUT:
-plot
(c) 2015, F Toso, R Garner, Centre for Future Air-Space Transportation
Technology, Univeristy of Strathclyde
%}

c_opt=param.LB+c0.*(param.UB-param.LB);
% optional scaling routine for engine and surfaces
[phases] = ExtractScalingPrameters (param, c_opt, phases);

loadconstants
[res] = extract_control_law (param, c_opt, phases);                          % initialize control matrix, time vector and x0 vector
[res] = propagate_trajectory (param, phases, const, res);

accelerations=[];    thermal=[];  forces =[];   controls =[]; Mn =[]; atmo=[];
index=0;
for ip = 1 : param.np                                                       % select phase
    for ine = 1 : phases(ip).ne                                             % select element
        index = index+1;                                                    % propagate trajectory
        for iicc=1:size(res(index).x,1)      
            [~, con]= phases(ip).dynamics(res(index).t(iicc), res(index).x(iicc,:), res(index).control, const, phases(ip));
            accelerations=[accelerations; con.acc]; %[accx, accz]
            thermal=[thermal; con.temp];            %[temperature, heat]
            forces=[forces; con.forces];            %[CL,CD,FT,q]
            controls=[controls; con.controls];
            Mn=[Mn; con.Mn];
           	atmo = [atmo;con.atmo];
        end
    end
end

plotroutines(phases, res, controls) % plot it
drawnow; pause(1);

advanced_plots (phases, res, forces, accelerations, Mn, thermal) % plot L/D, T, ACC, temp (option inside for partial selection of graph)
drawnow; pause(1);