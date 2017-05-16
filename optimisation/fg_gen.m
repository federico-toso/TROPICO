function [c0_opt] = fg_gen (c0, LB, UB, param, phases)
%{ 
[c0_opt] = fg_gen (c0, LB, UB, param, phases)

This function takes the first guess, uses the user defined algorithm to
find a better solution
 
INPUT
c0 = initial guess / generated vector for optimisation
LB = lower bound
UB = upper bound
param = general parameters structure
phases = different phases' custom parameters

OUTPUT
c0_opt = first guess refined

(c) 2015, F Toso, Centre for Future Air-Space Transportation Technology, Univeristy of Strathclyde
%}  

%% initial setup
c0_opt=c0;                              % optimizer vector is temporarily the first guess
usefmincon = 0; usega = 0;
if strcmp(param.fg{1},'lhs_fmincon')    % use matlab's fmincon and lhs to generate first guess
    usefmincon = 1;                 
elseif strcmp(param.fg{1},'ga')         % use matlab's ga to generate first guess
    usega = 1;
elseif strcmp(param.fg{1},'load')
    load(param.fg{2},'c0')
    c0_opt = c0;%(1:end-1);
    disp ('Previous solution loaded, Skipping first guess generation')
    disp('-------------------------------------------------------------------')
    return
elseif strcmp(param.fg{1},'load_noise')
%     load(param.fg{2},'c0')
    c0 = param.fg{2};
    options=optimset('Display', 'none' ,'MaxFunEvals',size(c0,2)*250,...
        'ScaleProblem','none','Algorithm', 'sqp','TolCon', 1e-3,'UseParallel',true);
    for inoise = -11:2:-1
        disp([datestr(now,'HH:MM dd-mmm'),' | Noise is 1e', num2str(inoise)])
        c0_noise = c0+10^(inoise)*c0.*(-0.5+rand(1,length(c0)));
        [c0,~,exitflag] = objconstr (c0_noise, LB, UB, options, param, phases);
        disp([datestr(now,'HH:MM dd-mmm'),' | New exitflag = ',num2str(exitflag)])
        [~,c,ceq] = compute_cost_const(c0, param, phases);
        if exitflag >=0 && max(c)<=1e-4 && max(abs(ceq))<=1e-4
            disp([datestr(now,'HH:MM dd-mmm'),' | Positive exitflag, minimum recovered, breaking'])
            c0_opt = c0;
            eval_plot (c0_opt, param, phases), drawnow, pause (0.5)
            disp ('Previous solution loaded, Skipping first guess generation')
            disp('-------------------------------------------------------------------')
            return
        end
    end
    usefmincon = 1;
    param.fg{2} = 1;
elseif strcmp(param.fg{1},'load_last')
    load(param.fg{2},'c0')
    disp ('Previous solution loaded, Skipping first guess generation')
    disp('-------------------------------------------------------------------')
    return
elseif strcmp(param.fg{1},'random')
    c0_opt = rand(size(c0));
    disp ('Random guess generated, Skipping first guess generation')
    disp('-------------------------------------------------------------------')
    return
elseif strcmp(param.fg{1},'lhs_feasibility_fmincon')
    original_cost = param.objfun;
    param.objfun = 'feasibility';
    param.fg{1} = {'lhs_fmincon'};
    param.keep_population=[];
    disp ('Starting feasibility anlysis')
    [c0_opt] = fg_gen (c0, LB, UB, param, phases);
    disp ('LHS_feasibility_fmincon first pass done')
    rmfield(param,'keep_population');
    population = c0_opt;
    param.objfun = original_cost;
    usefmincon = 1;
else                                    % skip first guess
    disp ('Skipping Direct multiple shooting first guess optimisation')
    disp('-------------------------------------------------------------------')
    return
end

disp ('Direct multiple shooting first guess optimisation started') 

if param.fg{2}~=1                       % timestep multiplier: >1 values reduces computation time     
    tstepm = param.fg{2};
else
    tstepm = 1;
end

% remove warnings
warning('off','all')
if param.parallel == 1
    pctRunOnAll  warning('off','all') %warning off for all the parallel workers
end

for ip=1:param.np % apply timestep multiplier for faster itegration
    phases(ip).tstep = phases(ip).tstep*tstepm;
%     param.odet = ceil(param.odet/tstepm);
end
disp (['First guess timestep: ',num2str(tstepm),'x'])

%% optimize control law

% initialize the first population, used for fmincon too
if ~exist('population','var') %i
    popsz = size(c0,2)+1;                                                       
    if isfield(param,'debug')
        popsz = 10; % use for quick debug or test
    end
    for ipo = 1:popsz                                                           % create the first population with LHS
        population = [ c0; bsxfun(@plus,bsxfun(@times,lhsdesign(popsz-1,size(c0,2)),(UB-LB)),LB)];
    end
else
    popsz = size(population,1);
end
exitflag = 0;  runs = 1;

% Matlab GA
if usega == 1                                                               % run matlab GA routine
    disp('Matlab GA started')
    timeval = tic;                                                          % start timing
    param.feq=0;                                                            % use inequality constraints
    gaoptions = gaoptimset('Display','iter','Generations',50,'PopulationSize',popsz,...
                'UseParallel',true,'TimeLimit',60*popsz,'InitialPopulation',population,...
                'Tolcon',1e-3);
    while exitflag<1 && runs <=3                                               % run GA algorithm *runs* times
        [c0_optGA,~,exitflag,~,population,scores] = ga_objconstr (c0, LB, UB, gaoptions, param, phases);   
        runs = runs+1;                                                      % run until solution find or max 3 times
    end
    if exitflag > 0                                                            % if GA succeded
        population = sortrows([population,scores],size(population,2)+1);    % reorder with scores
        population = population(:,1:end-1);                                 % remove scores
        [c0_opt] = find_best (c0_optGA, c0_opt, param, phases);    % check if solution is the new best
    else                                                                    % if ga failed, restart LHS population
        population = [ c0; bsxfun(@plus,bsxfun(@times,lhsdesign(popsz-1,size(c0,2)),(UB-LB)),LB)];
    end
    dureval = toc(timeval);                                                 % end timing
    disp(['Matlab GA completed, duration: ', num2str(dureval),' s'])
    disp('-------------------------------------------------------------------')
    population(1,:)=c0_opt;                                                 % first vector is now GA solution
    param.feq=1;
    exitflag=repmat(exitflag,popsz,1); % fix for refinement
end

% Fmincon
if usefmincon == 1
    disp([datestr(now,'HH:MM dd-mmm'),' | Matlab fmincon started on LHS population of ',num2str(popsz),' elements'])
    ttimeval = tic;                                                         % start timing
    options=optimset('Display', 'none' ,'MaxFunEvals',size(c0,2)*100,...
        'ScaleProblem','none','Algorithm', 'sqp','MaxSQPIter',5e3,'TolCon', 1e-4,'UseParallel',true);
    for ilhs = 1:popsz
        timeval = tic;                                                      % start partial timing, run optimisation
        [population(ilhs,:),~,exitflag(ilhs)] = objconstr (population(ilhs,:), LB, UB, options, param, phases);
        save('population','population','exitflag');
        dureval = toc(timeval);                                             % end partial timing
        completiontext=[num2str(ilhs/popsz*100,'%5.1f'),'%'];
        elapsed_time =toc(ttimeval);
        timeleft =  elapsed_time / ilhs*popsz - elapsed_time; % * (popsz/(popsz+4));
        disp([datestr(now,'HH:MM dd-mmm'),...
            ' | ',completiontext,' completed',...
            ' | exitflag: ',num2str(exitflag(ilhs),'%+2.0f'),...
            ' | duration: ',num2str(round(dureval)),' s',...
            ' | ETC: ',datestr(now+1/24/60/60*timeleft,'HH:MM dd-mmm')]);
    end
    disp(['The solution converged ',num2str(sum(exitflag>=0)),' times'])
    tdureval = toc(ttimeval);                                               % end total timing
    disp(['Matlab fmincon evaluation of LHS population completed, duration: ', num2str(round(tdureval/60)),' min'])
    disp('-------------------------------------------------------------------')
end


disp('Direct multiple shooting first guess optimisation completed')

%% refine solution for different timesteps
disp ('Refining solution for the correct timestep...')
for ip=1:param.np
    phases(ip).tstep = phases(ip).tstep/tstepm;
end
options=optimset('Display', 'none' ,'MaxFunEvals',size(c0,2)*250,...
        'ScaleProblem','none','Algorithm', 'sqp','TolCon', 1e-3,'UseParallel',true);
num_succ=0;
exitflag_new = exitflag; %uncomment if noise to everything
if param.fg{2}~=1
    for ifinalize = 1:popsz
%         if exitflag(ifinalize)>=0
            num_succ=num_succ+1;
                disp([datestr(now,'HH:MM dd-mmm'),' | Starting on element #',num2str(ifinalize),', success #',num2str(num_succ),' that had exitflag = ', num2str(exitflag(ifinalize))])
                [c_opt_final(num_succ,:),fval_final(num_succ),exitflag_new(num_succ)] = objconstr (population(ifinalize,:), LB, UB, options, param, phases);
                disp([datestr(now,'HH:MM dd-mmm'),' | New exitflag = ',num2str(exitflag_new(num_succ))])
                if exitflag_new(num_succ) < 0
                    disp([datestr(now,'HH:MM dd-mmm'),' | Starting tentative recovery with noise'])
                    for inoise = -7:2:-3
                        disp([datestr(now,'HH:MM dd-mmm'),' | Noise is 1e', num2str(inoise)])
                        c0_noise = population(ifinalize,:)+10^(inoise)*population(ifinalize,:).*(-0.5+rand(1,length(c0)));
                        [c_opt_final(num_succ,:),fval_final(num_succ),exitflag_new(num_succ)] = objconstr (c0_noise, LB, UB, options, param, phases);
                        disp([datestr(now,'HH:MM dd-mmm'),' | New exitflag = ',num2str(exitflag_new(num_succ))])
                        if exitflag_new(num_succ) >=0
                            disp([datestr(now,'HH:MM dd-mmm'),' | Positive exitflag, minimum recovered, breaking'])
                            break
                        end
                    end
                end
%         end
    end
else
    disp('Timestep multiplier is 1x, skipping refinement...')
    c_opt_final=population;
    num_succ=popsz;
end

% find best value
if ~isfield(param,'keep_population')
    disp('Finding best solution among converged...')
    for ifinal=1:num_succ
        if exitflag_new(ifinal) >= 0
            disp(['Comparing temporary best solution with case ',num2str(ifinal)])
            [c0_opt] = find_best (c_opt_final(ifinal,:), c0_opt, param, phases);
            if c_opt_final(ifinal,:)==c0_opt
                eval_plot (c0_opt, param, phases)
                drawnow; pause(0.05);
            end
        end
    end
end

cd('./results');
save(['firstguess_', datestr(now, 'yyyymmddTHHMMSS')])
cd('../');

cd('./results');                                                  % save the first guess results
save('firstguess')
cd('../');
disp ('Direct multiple shooting first guess generated')
disp('-------------------------------------------------------------------')

%% final checks
if size(c0_opt,2) ~= size(c0,2)                                             % check if there is a size error
    error('size of the output of fg_gen.m is different from the input')
end

for chep = 1:size(c0,2)                                                     % check if all values are between bounds
    if c0_opt(chep)>UB(chep) || c0_opt(chep)<LB(chep)
        disp(['ne: ',num2str(ne),' nc: ',num2str(nc)])
        disp(['index: ',num2str(chep),' c0_opt(chep)= ',num2str(c0_opt(chep)),...
            ' LB(chep)= ',num2str(LB(chep)),' UB(chep)= ',num2str(UB(chep)),])
        warning('some values of the output of fg_gen.m are out of bounds')
    end
end

end



%% function to evaluate if the new optimised solution improved
function [c0_opt] = find_best (new_c0, old_c0, param, phases)
    [old_cost_function,oc,oceq] = compute_cost_const(old_c0, param, phases);
    if isempty(oc); oc=0; end;                              % fill constraints if empty
    if isempty(oceq); oceq=0; end;                          % fill constraints if empty
    if max(oc)>=1e-2 || max(abs(oceq))>=1e-2                % check if previous solution was good
        old_cost_function = 1e6;                            % if not, put placeholder cost function
    end                                                     % evaluate the new optimisation vector
    [new_cost_function,c,ceq] = compute_cost_const(new_c0, param, phases);
    if isempty(c); c=0; end;                                % fill constraints if empty
    if isempty(ceq); ceq=0; end;                            % fill constraints if empty
    if max(c)<=1e-2 && max(abs(ceq))<=1e-2 && new_cost_function<=old_cost_function
        c0_opt = new_c0;                                    % if new best solution found, save it
        eval_plot (c0_opt, param, phases)       % plot result
        cd('./results');
        save('firstguess_temp')                             % save temporary best solution found
        cd('../');
        disp('New best first guess found!  <<<')
    else
        c0_opt = old_c0;                                    % if not, keep old solution
        if max(c)>1e-2 || max(abs(ceq))>1e-2
            disp('false positive exitflag found, discarding')
%             if max(c)>1e-2
%                 disp('c:')
%                 c
%             end
%             if max(abs(ceq))>1e-2
%                 disp('ceq:')
%                 ceq
%             end
        end
    end
    drawnow; pause(0.5);
end