function run_TROPICO(param, phases)
%{
///////////////////////////////// to be updated/////////////////////
main_optimisation(InitialParameters,AscTargetParameters,ref)

Main function to generate an optimal control law based on the angle of
attack and mission time for the rocket segment of an ascent trajectory,
configured for the Skylon C1 vehicle based on OPTRAJ.
The optimisation varies the control law with an aim to maximise the final
vehicle mass into orbit. The optimisation is based on the built-in MATLAB
function 'fmincon' using a SQP algorithm with an additional multiple
shooting approach.

INPUT
lat0 = Initial latitude at the start of the rocket phase (rad N)
lon0 = Initial longitude at the start of the rocket phase (rad E)
horb = Target final altitude (m)
t_e = Initial guess for ascent phase duration (s)

FUNCTIONS CALLED
2- Optimisation cost function, optrajCostFunction_MS.m
3- Optimisation constraints function, optrajConstraints_MS.m
4- Ordinary differential equations of motion for modelling the ascent trajectory, StrathEOM.m
///////////////////////////////// to be updated/////////////////////

(c) 2012, C Maddock & E Minisci, Centre for Future Air-Space Transportation
 Technology, Univeristy of Strathclyde
%}

%% initial post
if isfield(param,'debug')
    disp('TROPICO is running in debug mode, plotting and first guess complexity reduced')
end

disp('Ascent Optimization started')

%% Initialise some variables
loadconstants                           % Load constants (remove in the future)
tottime=tic;                            % start timer
if exist('results','dir')~=7
    mkdir results;
end

param.ns=numel(phases(1).x0);            % Number of elements in the state vector, i.e, [1-h, 2-v, 3-gamma, 4-chi, 5-lat, 6-lon, 7-m]                                                   
param.nv=size(phases(1).cbounds,1);      % Number of elements in the control vector [pitch, throttle, bank]
param.np =numel(phases);
input_updater;

%% Parallel pool settings
poolobj = gcp('nocreate');
if param.parallel > 0
    if isempty(poolobj)
        if param.parallel == 1
            parpool;
        else
            parpool(param.parallel);
        end
    end
elseif param.parallel == 0
    delete(poolobj);
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
end

%% Transcription of opt vector and scaling between 0 & 1 to better tune the optimiser

disp('Encoding bounds...')
[LB, UB, param, phases] = generate_bounds (param, phases);

%% generation of the first guess / population
disp('Generating initial guess/population...')
[c0_population, param, phases] = generate_guess (param, phases);
%%%%%
% c0_population=c0_population*0+1;
%%%%
cd('./results');
save('initial_population.mat')
cd('../');

%% analysis of the population to select the best candidate
if size(c0_population,1)~=1
    disp('Starting individual/population analysis...')
    fgtime = tic();
    [c0, ranked_population] = run_optimisation (c0_population, param, phases);
    fg_time = toc(fgtime);
    disp(['First guess search duration: ', num2str(round(fg_time/60)),' min'])
else
    c0 = c0_population;
    ranked_population = c0;
end

cd('./results');
save('optimised_population.mat')
delete('initial_population.mat')
cd('../');

%% final optimization

% Optimisation settings:
options=optimset('Algorithm', 'sqp', ...                         % declare algorithm to use, here set to sequential quadratic programming (SQP)
    'Display', 'iter-detailed', ...                                         %  displays output at each iteration, and gives the technical exit message.
    'TolCon', 1e-4, ...                                                     % tolerance on the constraint violation %'TolX', 1e-9, ...         % tolerance on the design vector     %
    'MaxSQPIter',5e3, ...                                                      % maximum number of SQP iterations allowed
    'MaxFunEvals',size(c0,2)*1000, ...                                      % maximum number of function evaluations allowed
    'ScaleProblem','none', ...                                    %  causes the algorithm to normalize all constraints and the objective function    
    'UseParallel',true...
    );  %     'TolFun', 1e-6, ...                                                     % termination tolerance on the function value


disp('Direct multiple shooting fmincon optimisation started')

phases_refined = phases;

% if size(c0_population,1)~=1
%     for ip = 1:param.np
%         phases_refined(ip).ode=@ode5;
%         phases_refined(ip).tstep=phases(ip).tstep/2;
%     end
% end
% [c0, ~,exitflag,~,~,~,~] = objconstr (c0, LB, UB, options, param, phases_refined);

%% Integrate trajectory based on the inital, final times, 
% initial state vector and optimised control law
% subject to conditions given in events_alt.m, i.e., final altitude         

c0_unscaled=param.LB+c0.*(param.UB-param.LB);
[res] = extract_control_law (param, c0_unscaled, phases_refined);
% scaling
[phases_refined] = ExtractScalingPrameters (param, c0_unscaled, phases_refined);
% integration
[res] = propagate_trajectory (param, phases_refined, const, res);
% computation of cost and constraints
[cost_function,c,ceq] = compute_cost_const(c0, param, phases);

accelerations=[];    thermal=[];  forces =[];   controls =[]; Mn =[]; atmo=[];
index=0;
for ip = 1 : param.np                                                       % select phase
    for ine = 1 : phases_refined(ip).ne                                             % select element
        index = index+1;                                                    % propagate trajectory
        for iicc=1:size(res(index).x,1)      
            [~, con]= phases_refined(ip).dynamics(res(index).t(iicc), res(index).x(iicc,:), res(index).control, const, phases_refined(ip));
            accelerations=[accelerations; con.acc]; %[accx, accz]
            thermal=[thermal; con.temp];            %[temperature, heat]
            forces=[forces; con.forces];            %[CL,CD,FT,q]
            controls=[controls; con.controls];
            Mn=[Mn; con.Mn];
           	atmo = [atmo;con.atmo];
        end
    end
end

% end timer
totaltime=toc(tottime);
disp(['Total duration: ', num2str(round(totaltime/60)),' min'])

%% save variables to a mat file
if exist('exitflag','var')
    if exitflag>0
        nm = ['success_', datestr(now, 'yyyymmddTHHMMSS')];
    else
        nm = ['fail_', datestr(now, 'yyyymmddTHHMMSS')];
    end
else
    nm = ['case_', datestr(now, 'yyyymmddTHHMMSS')];
end
cd('./results');
save(nm)
delete('optimised_population.mat')
cd('../');
disp('Optimization output saved')

%% Plotting Routines

plotroutines(phases_refined, res, controls)
advanced_plots(phases_refined, res, forces, accelerations, Mn, thermal)
if ~isfield(param,'debug')
    finalplot(param, phases_refined, res)
end
drawnow; pause(0.5);
    
%% Display MDO statistics 

format short g
if param.ScalingActive == 1
    disp('-------------------------------------------------------------------')
    k=length(c0);
    if isfield(param,'Scale_S2_Surf')
        if param.Scale_S2_Surf(1)-param.Scale_S2_Surf(2)~=0
            disp(['Stage 2 aerodynamic surface scaling = ',num2str(c0_unscaled(k))])
            k = k-1;
        end
    end
    if isfield(param,'Scale_S1_Surf')
        if param.Scale_S1_Surf(1)-param.Scale_S1_Surf(2)~=0
            disp(['Stage 1 aerodynamic surface scaling = ',num2str(c0_unscaled(k))])
            k = k-1;
        end
    end
    if isfield(param,'Scale_S2_Engine')
        if param.Scale_S2_Engine(1)-param.Scale_S2_Engine(2)~=0
            disp(['Stage 2 Engine scaling = ',num2str(c0_unscaled(k))])
            k = k-1;
        end
    end
    if isfield(param,'Scale_S1_Engine')
        if param.Scale_S1_Engine(1)-param.Scale_S1_Engine(2)~=0
            disp(['Stage 1 Engine scaling = ',num2str(c0_unscaled(k))])
        end
    end    
end  
if isfield(param,'Position_S2_GTOW')
    disp(['Stage 2 GTOW = ',num2str(c0_unscaled(param.Position_S2_GTOW)),' kg'])
end

%% Display trajectory statistics

disp('-------------------------------------------------------------------')
disp('Final Conditions')
disp(['Altitude = ',num2str(res(end).x(end,1)),' [m]'])
disp(['Velocity = ',num2str(res(end).x(end,2)),' [m/s]'])
disp(['Longitude = ', num2str(res(end).x(end,6)*180/pi),' [deg]'])
disp(['Latitude = ', num2str(res(end).x(end,5)*180/pi),' [deg]'])
disp(['Flightpath angle = ',num2str(res(end).x(end,3)*180/pi),' [deg]'])
disp(['Heading = ',num2str(res(end).x(end,4)*180/pi),' [deg]'])
disp(['Mass = ',num2str(res(end).x(end,7)),' [kg]'])
disp(['Mission duration = ',num2str(res(end).t(end)),' [s]'])
disp('-------------------------------------------------------------------')
