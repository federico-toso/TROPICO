% update input file to the last version

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
phases(end).odet = [];
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
phases(end).controls_no_match = [];
phases(end).nu = [];
phases(end).continue_from = [];
phases(end).ceq=[];

phases = phases(1:param.np);