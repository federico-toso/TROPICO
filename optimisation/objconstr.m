function [c_opt, fval,exitflag,output,lambda,grad,hessian] = objconstr...
    (c0, LB, UB, options, param, phases)

xLast = [];
myf = [];
myc = [];
myceq = [];

fun = @objfun;
cfun = @constr;

% rescale relative constraints
[~,c,ceq] = compute_cost_const(c0, param, phases);
% scale_c = abs(c);
% scale_c(scale_c<=1)=1;
% scale_ceq = ceq;
% scale_ceq(scale_ceq<=1)=1;
options=optimset(options, 'TolCon', options.TolCon/max([max(c(c>0)),max(abs(ceq)),1]));

% run optimisation
[c_opt, fval,exitflag,output,lambda,grad,hessian] = fmincon(...
    fun, c0, [],[],[],[], LB, UB, ...
    cfun, options, param, phases);

% check absolute constraints <1e-3
[exitflag] = absolute_exitflag (c_opt, param, phases, exitflag);
	
	function obj = objfun(c_opt, param, phases)
        if ~isequal(c_opt,xLast) % Check if computation is necessary
			[myf,myc,myceq] = compute_cost_const(c_opt, param, phases);
			xLast = c_opt;
        end
        obj = myf;
    end		
		
	function [c, ceq] = constr(c_opt, param, phases)
		if ~isequal(c_opt,xLast) % Check if computation is necessary
			[myf, myc, myceq] = compute_cost_const(c_opt, param, phases);
			xLast = c_opt;
		end
		c = myc;
		ceq = myceq;
    end

    function [exitflag] = absolute_exitflag (c_opt, param, phases,exitflag)
        [c, ceq] = constr(c_opt, param, phases);

        max_c=max(c) ;
        max_ceq=max(abs(ceq)) ;
        if isempty(max_c)
            max_c=0;
        end
        if isempty(max_ceq)
            max_ceq = 0;
        end
        
        if max_c < 1e-4 && max_ceq < 1e-4 && exitflag<1
            exitflag = 6;
        end
    end

end