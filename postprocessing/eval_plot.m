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

plotroutines(phases, res) % plot it
drawnow; pause(1);

advanced_plots(phases, res) % plot L/D, T, ACC, temp (option inside for partial selection of graph)
drawnow; pause(1);