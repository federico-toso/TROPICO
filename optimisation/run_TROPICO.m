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

%% display initial log
% disp('-------------------------------------------------------------------')
% disp('Initial Conditions: ')
% disp(['Altitude = ',num2str(initial.altitude),' [m]'])
% disp(['Velocity = ',num2str(initial.velocity),' [m/s]'])
% disp(['Longitude = ',num2str(rad2deg(initial.longitude)),' [deg]'])
% disp(['Latitude = ',num2str(rad2deg(initial.latitude)),' [deg]'])
% disp(['Flightpath angle = ',num2str(rad2deg(initial.flightPathAngle)),' [deg]'])
% disp(['Heading = ',num2str(rad2deg(initial.headingAngle)),' [deg]'])
% disp(['Mass = ',num2str(initial.mass),' [kg]'])
% disp('-------------------------------------------------------------------')
% disp('Target Conditions: ')
% disp(['Altitude = ',num2str(target.altitude),' [m]'])
% disp(['Velocity = ',num2str(target.velocity),' [m/s]'])
% disp('-------------------------------------------------------------------')

% disp('Simulation setup: ')
% disp(['Multiple shooting elements = ',num2str(param.ne)])
% disp(['Control points for each phase = ',num2str(param.nc)])
% disp(['Control vector length = ',num2str(param.ne*(1+3*param.nc)+(param.ne-1)*(7)),' elements'])
% disp(['Simulation timestep = ',num2str(param.tstep),' [s]'])
% disp('-------------------------------------------------------------------') 
if isfield(param,'debug')
    disp('TROPICO is running in debug mode, plotting and first guess complexity reduced')
end

disp('Ascent Optimization started')

%% Initialise some variables
loadconstants                           % Load constants (remove in the future)
tottime=tic;                                     % start timer
if exist('results','dir')~=7
    mkdir results;
end

param.ns=7;                             % Number of elements in the state vector, i.e, [1-h, 2-v, 3-gamma, 4-chi, 5-lat, 6-lon, 7-m]                                                   
param.nv=3;                             % Number of elements in the control vector [pitch, throttle, bank]
param.feq = 1;                          % parameter used to choose beteween equality or inequality constraints

poolobj = gcp('nocreate');
if param.parallel == 1
    if isempty(poolobj)
        parpool;
    end
elseif param.parallel == 0
    delete(poolobj);
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
end


%% To better tune to optimiser, the values of the state vector are scaled to be closer to unity (i.e., [0, 1]) with minimum & maximum values declared
LB=[];
UB=[];
c0=[];
param.np = numel(phases);
input_updater;                                                              % older input files get updated
if ~isfield(param,'bias_obj_const')
    param.bias_obj_const = [1,1];
end

for ip = 1 : param.np                                                       % select phase
    phases(ip).nu = size(phases(ip).cbounds,1);                             % Number of elements in the control vector, e.g., [angle of attack, throttle, bank]
    for ine = 1 : phases(ip).ne                                             % select element
%         add control values, and bounds
        for iv = 1:phases(ip).nu
            if phases(ip).cbounds(iv,1)~=phases(ip).cbounds(iv,2)
                c0 = [c0, phases(ip).fg(iv)         *ones(1,phases(ip).nc)];    % first guess
                LB = [LB, phases(ip).cbounds(iv,1)  *ones(1,phases(ip).nc)];    % lower bound
                UB = [UB, phases(ip).cbounds(iv,2)  *ones(1,phases(ip).nc)];    % upper bound
            end
        end
%         add time variable
        phases(ip).tof(1) = max(phases(ip).tof(1),1);
        c0 = [c0, (phases(ip).tof(1)+phases(ip).tof(2))/2/phases(ip).ne];   %first guess is average between bounds
        LB = [LB, phases(ip).tof(1)/phases(ip).ne];
        UB = [UB, phases(ip).tof(2)/phases(ip).ne];
%         add state variables
        if isempty(phases(ip).x0_no_opt) || ine ~=1                         % if x0_no_opt is not defined, add all
            LB = [LB, phases(ip).xbounds(:, 1)'];
            UB = [UB, phases(ip).xbounds(:, 2)'];
            c0 = [c0, (phases(ip).xbounds(:, 2)+phases(ip).xbounds(:, 1))'/2];
        else                                                                % if not, add the specified ones
            for ins = 1 : size(phases(ip).xbounds,1)                        % for each state variable
                if phases(ip).x0_no_opt(ins)~=1 
                    c0 = [c0, phases(ip).x0(ins)];                          % if it's not fixed, add value
                    LB = [LB, phases(ip).xbounds(ins, 1)];
                    UB = [UB, phases(ip).xbounds(ins, 2)];
                end
            end
        end
        if ip>1 && ine ==1                                                  % in case of staging, don't optimise mass, remove last element
            if phases(ip).vehicle.gtow~=phases(ip-1).vehicle.gtow
                c0=c0(1:end-1);    LB=LB(1:end-1);    UB=UB(1:end-1);
            	if isfield(param,'Scale_S2_GTOW') && ~exist('gtow_s2_opt','var')   % second stage is optimized, add values for it
                    if param.Scale_S2_GTOW(1)-param.Scale_S2_GTOW(2)~=0             % if optimisation bounds are different
                        gtow_s2_opt=1;
                        c0 = [c0, mean(param.Scale_S2_GTOW)];
                        LB = [LB, param.Scale_S2_GTOW(1)];
                        UB = [UB, param.Scale_S2_GTOW(2)];
                        param.Position_S2_GTOW = length(c0);
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
%% additional extra optimisation variables for scaling of forces (L,D, T, mp)
if isfield(param,'Scale_S1_Engine') || isfield(param,'Scale_S2_Engine') || isfield(param,'Scale_S1_Surf') ||...
        isfield(param,'Scale_S2_Surf') || isfield(param,'Scale_S2_GTOW')
    param.ScalingActive=1;
    [c0, LB, UB, phases] = AddScalingPrameters (c0, LB, UB, param, phases);
else
    param.ScalingActive=0;
end

%% new transcription
param.LB=LB;
param.UB=UB;
c0=(c0-LB)./(UB-LB);
UB=ones(1,numel(UB));
LB=0*UB;


% opt vector checks
if sum(param.LB==param.UB)>0
    error(['There are ',num2str(sum(x_scale==0)),' variables with the same upper and lower bounds.'])
end
if sum(c0<=LB)>0 || sum(c0>=UB)>0
    c0(c0<=LB)=LB(c0<=LB)+1e-9;
    c0(c0>=UB)=UB(c0>=UB)-1e-9;
    disp('First guess is out of bounds, values shifted')
end


%% generation fo the first guess
fgtime = tic();
[c0fg] = fg_gen (c0, LB, UB, param, phases);
fg_time = toc(fgtime);
disp(['First guess search duration: ', num2str(round(fg_time/60)),' min'])
param.bias_obj_const = [1,1];
%% Run optimisation
% Optimisation settings:
options=optimset('Algorithm', 'sqp', ...                         % declare algorithm to use, here set to sequential quadratic programming (SQP)
    'Display', 'iter-detailed', ...                                         %  displays output at each iteration, and gives the technical exit message.
    'TolCon', 1e-4, ...                                                     % tolerance on the constraint violation %'TolX', 1e-9, ...         % tolerance on the design vector     %
    'MaxSQPIter',5e3, ...                                                      % maximum number of SQP iterations allowed
    'MaxFunEvals',size(c0,2)*1000, ...                                      % maximum number of function evaluations allowed
    'ScaleProblem','none', ...                                    %  causes the algorithm to normalize all constraints and the objective function    
    'UseParallel',true...
    );  %     'TolFun', 1e-6, ...                                                     % termination tolerance on the function value

% final optimization 
disp('Direct multiple shooting fmincon optimisation started')
[c0, ~,exitflag,~,~,~,~] = objconstr (c0fg,LB, UB, options, param, phases);

%% Integrate trajectory based on the inital, final times, initial state vector and optimised control law
% subject to conditions given in events_alt.m, i.e., final altitude

c0_unscaled=param.LB+c0.*(param.UB-param.LB);
[res] = extract_control_law (param, c0_unscaled, phases);
% scaling
[phases] = ExtractScalingPrameters (param, c0_unscaled, phases);
% integration
[res] = propagate_trajectory (param, phases, const, res);

totaltime=toc(tottime);                                                                         % end timer
disp(['Total duration: ', num2str(round(totaltime/60)),' min'])
%% save variables to a mat file ------------------------------------------
if exitflag>0
    nm = ['success_', datestr(now, 'yyyymmddTHHMMSS')];
else
    nm = ['fail_', datestr(now, 'yyyymmddTHHMMSS')];
end
cd('./results');
save(nm)
disp('Optimization output saved')
cd('../');

%% Plotting Routines -----------------------------------------------------
plotroutines(phases, res)
advanced_plots(phases, res)
if ~isfield(param,'debug')
    finalplot(param, phases, res)
end
drawnow; pause(0.5);
    
%% Display MDO statistics -------------------------------------------------
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
%% Display trajectory statistics ------------------------------------------
% disp('-------------------------------------------------------------------')
% disp('Final Conditions')
% disp(['Altitude = ',num2str(res(end).x(end,1)),' [m]'])
% disp(['Velocity = ',num2str(res(end).x(end,2)),' [m/s]'])
% disp(['Longitude = ', num2str(res(end).x(end,6)*180/pi),' [deg]'])
% disp(['Latitude = ', num2str(res(end).x(end,5)*180/pi),' [deg]'])
% disp(['Flightpath angle = ',num2str(res(end).x(end,3)*180/pi),' [deg]'])
% disp(['Heading = ',num2str(res(end).x(end,4)*180/pi),' [deg]'])
% disp(['Mass = ',num2str(res(end).x(end,7)),' [kg]'])
% disp(['Mission duration = ',num2str(res(end).t(end)),' [s]'])
% disp('-------------------------------------------------------------------')