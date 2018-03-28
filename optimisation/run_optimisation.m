function [c0_opt, ranked_population] = run_optimisation  (c0_population, param, phases)

popsz=size(c0_population,1);
UB=ones(1,size(c0_population,2));
LB=0*UB;
exitflag = -1*ones(popsz,1);

disp([datestr(now,'HH:MM dd-mmm'),' | Matlab fmincon started on population of ',...
    num2str(popsz),' individuals'])
ttimeval = tic;                                                         % start timing
options=optimset('Display', 'none' ,'MaxFunEvals',size(c0_population,2)*100,...
    'ScaleProblem','none','Algorithm', 'sqp','MaxSQPIter',5e3,'TolCon', 1e-4,'UseParallel',true);
for i = 1:popsz
    timeval = tic;                                                      % start partial timing, run optimisation
    [c0_population(i,:),~,exitflag(i)] = objconstr (c0_population(i,:), LB, UB, options, param, phases);
    save('population','c0_population','exitflag');
    dureval = toc(timeval);                                             % end partial timing
    completiontext=[num2str(i/popsz*100,'%5.1f'),'%'];
    elapsed_time =toc(ttimeval);
    timeleft =  elapsed_time / i*popsz - elapsed_time; % * (popsz/(popsz+4));
    disp([datestr(now,'HH:MM dd-mmm'),...
        ' | ',completiontext,' completed',...
        ' | exitflag: ',num2str(exitflag(i),'%+2.0f'),...
        ' | duration: ',num2str(round(dureval)),' s',...
        ' | ETC: ',datestr(now+1/24/60/60*timeleft,'HH:MM dd-mmm')]);
    cd('./results');
    save('population_temp.mat')
    cd('../');
end
disp(['The solution converged ',num2str(sum(exitflag>=0)),' times'])
tdureval = toc(ttimeval);                                               % end total timing
disp(['Matlab fmincon evaluation of LHS population completed, duration: ', num2str(round(tdureval/60)),' min'])
disp('-------------------------------------------------------------------')

%% find best value
disp('Finding best solution among converged...')
c0_opt=LB;
successes = 0;
for i=1:popsz
    if exitflag(i) >= 0
        disp(['Checking individual # ',num2str(i)])
        [c0_opt] = find_best (c0_population(i,:), c0_opt, param, phases);
        if isequal(c0_population(i,:),c0_opt)
            successes = successes+1;
            if ~isfield(param,'debug')
                eval_plot (c0_opt, param, phases)
                drawnow; pause(0.05);
            end
        end
    end
end
%% if no individual satisfies the constraints, pick lowest constraint violation
if successes == 0
    disp('No solutions converged, choosing the case with lowest constraint violation...')
    max_violation = 1e3;
    for i=1:popsz
        [~,c,ceq] = compute_cost_const(c0_population(i,:), param, phases);
        new_max_violation = max(max(c),max(abs(ceq)));
        if new_max_violation<max_violation
            disp('New best first guess found!  <<<')
            c0_opt = c0_population(i,:);
            max_violation = new_max_violation;
        end
    end
end

%% TO ADD: ORDER POPULATION, right now this just stores all optimised cases
ranked_population = c0_population;

cd('./results');
delete('population_temp.mat')
cd('../');

end

%% function to evaluate if the new optimised solution improved
function [c0_opt] = find_best (new_c0, old_c0, param, phases)
    [old_cost_function,oc,oceq] = compute_cost_const(old_c0, param, phases);
    if isempty(oc); oc=0; end;                              % fill constraints if empty
    if isempty(oceq); oceq=0; end;                          % fill constraints if empty
    if max(oc)>=1e-4 || max(abs(oceq))>=1e-4                % check if previous solution was good
        old_cost_function = 1e6;                            % if not, put placeholder cost function
    end                                                     % evaluate the new optimisation vector
    [new_cost_function,c,ceq] = compute_cost_const(new_c0, param, phases);
    if isempty(c); c=0; end;                                % fill constraints if empty
    if isempty(ceq); ceq=0; end;                            % fill constraints if empty
    if max(c)<=1e-4 && max(abs(ceq))<=1e-4 && new_cost_function<=old_cost_function
        c0_opt = new_c0;                                    % if new best solution found, save it
        if ~isfield(param,'debug')
            eval_plot (c0_opt, param, phases)       % plot result
        end
        disp('New best first guess found!  <<<')
    else
        c0_opt = old_c0;                                    % if not, keep old solution
        if max(c)>1e-4 || max(abs(ceq))>1e-4
            disp('false positive exitflag found, discarding')
        end
    end
    drawnow; pause(0.5);
end
