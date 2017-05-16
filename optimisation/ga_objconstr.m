function [c_opt, fval,exitflag,output,population,scores] = ga_objconstr...
    (c0, LB, UB, options, param, phases)

xLast = [];
myf = [];
myc = [];
myceq = [];
fun = @(x)objfun(x, param, phases);
cfun = @(x)constr(x, param, phases);

[c_opt, fval,exitflag,output,population,scores] = ga(...
    fun, size(c0,2), [],[],[],[], LB, UB, ...
    cfun, options);

    function obj = objfun(c_opt, param, phases)
        if ~isequal(c_opt,xLast) % Check if computation is necessary
            [myf,myc,myceq] = compute_cost_const(c_opt, param, phases);
            xLast = c_opt;
        end
        obj = myf;
    end


    function [c,ceq] = constr(c_opt, param, phases)
        if ~isequal(c_opt,xLast) % Check if computation is necessary
            [myf, myc, myceq] = compute_cost_const(c_opt, param, phases);
            xLast = c_opt;
        end
        c = myc;
        ceq =myceq;
    end
end