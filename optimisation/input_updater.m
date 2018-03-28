%% update missing phase data with previous phase info
updatelist = {'vehicle','ne','nc','xbounds','cbounds','tof',...
    'ode','tstep','atmo','dynamics','aero','prop','thermal','inte','dist'};
for ip = (2:param.np)
    for ipr = 1:size(updatelist,1)
        checkstr = ['phases(',num2str(ip),').',cell2mat(updatelist(ipr))];
        if isempty(eval(checkstr))
            prev_value = ['phases(',num2str(ip-1),').',cell2mat(updatelist(ipr))];
            eval([checkstr,' = ',prev_value]);
        end
    end
end

%% create empty fields if not used
phases(end+1).vehicle.gtow = [];
phases(end).vehicle.m0 = [];
phases(end).vehicle.Sgross = [];
phases(end).vehicle.Rnose = [];
phases(end).vehicle.Rwing = [];
phases(end).vehicle.Epsilon = [];
phases(end).vehicle.Tangle = [];
phases(end).ne = [];
phases(end).nc = [];
phases(end).x0 = [];
phases(end).xbounds = [];
phases(end).cbounds = [];
phases(end).fg = [];
phases(end).tof = [];
phases(end).ode = [];
phases(end).testep = [];
phases(end).atmo = [];
phases(end).dynamics = [];
phases(end).prop = [];
phases(end).aero = [];
phases(end).thermal = [];
phases(end).inte = [];
phases(end).dist = [];
phases(end).xf = [];
phases(end).xf_fixed = [];
phases(end).xf_range = [];
phases(end).controls_no_match = [];
phases(end).nu = [];
phases(end).continue_from = [];
phases(end).ceq=[];
phases(end).control_tol=[];
phases(end).state_tol=[];

phases = phases(1:param.np);