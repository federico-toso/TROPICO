%{
Main function to optimize the mission, preset to optimise a trajectory

(c) 2015, F Toso, R Garner, Centre for Future Air-Space Transportation
Technology, Univeristy of Strathclyde
%}
%% initialize workspace --------------------------------------------------
% clc
clearvars
close all

%% Set initial conditions-------------------------------------------------
addpath(genpath(cd))
loadconstants

% Defintion of phases
% ------------------------------------------------------------------------
% WARNING: Phases MUST be defined by grouping the same stage together
% ------------------------------------------------------------------------
i=0;

%% Rocket ascent s1
i = i+1;
phases(i).vehicle = load_vehicle('falcon9_S1+2');
phases(i).atmo = @atmo_ext_ISA;
phases(i).dynamics = @dyn_ECEF2;
phases(i).prop = @prop_rocket_300s_1MN;
phases(i).aero = @aero_simplerocket;
phases(i).thermal = @temp_convradbalance;

phases(i).ne   = 2;             % number of elements per phase (for multi-shooting)
phases(i).nc   = 7;            % number of control nodes for each element within the phase
phases(i).inte = @interp_pch_mex;
phases(i).dist = @linspace; %@dist_chebyshev_mod;

phases(i).x0 =         [1000,  100, 85*pi/180, pi/2, 28.6083*pi/180, -80.6041*pi/180, 550e3-25*(7.6e6/282/9.8065)];
phases(i).x0_fixed =   [1, 1, 0, 0, 1, 1, 1];                              % id of fixed starting values

phases(i).xbounds = [10, 100e3;...            % altitude
                     10, 8e3; ...           % velocity
                     -pi/2, pi/2; 0, 2*pi; ...                                  % FPA, heading
                     -pi/2, pi/2; -pi, pi; ...       % lat, lon
                     phases(i).vehicle.m0, phases(i).vehicle.gtow+1];                           % mass
phases(i).cbounds = [-20*pi/180, 20*pi/180; 0.75, 1; 0, 0]; % [min, max] for alpha, throttle, beta

phases(i).tof = [60; 240];          % bounds for time of flight for the phase
phases(i).ode = @ode4;
phases(i).tstep= 1;                 % minimum timestep [seconds]




% Rocket ascent s2
i = i+1;
phases(i).vehicle = load_vehicle('falcon9_S2');
phases(i).atmo = @atmo_ext_ISA;
phases(i).dynamics = @dyn_ECEF2;
phases(i).prop = @prop_rocket_350s_1MN;
phases(i).aero = @aero_simplerocket;
phases(i).thermal = @temp_convradbalance;

phases(i).ne   = 2;             % number of elements per phase (for multi-shooting)
phases(i).nc   = 7;            % number of control nodes for each element within the phase
phases(i).inte = @interp_pch_mex;
phases(i).dist = @linspace; %@dist_chebyshev_mod;

phases(i).xf =          [200e3, 2e3, 0, pi/2, 0, 0, 50e3];
phases(i).xf_fixed =    [1, 0, 0, 0, 0, 0, 0];                              % id of target  value
phases(i).ceq = {'ceq_circ',0};

phases(i).xbounds = [30e3, 200e3;...            % altitude
                     1e3, 8e3; ...           % velocity
                     -pi/2, pi/2; 0, 2*pi; ...                                  % FPA, heading
                     -pi/2, pi/2; -pi, pi; ...       % lat, lon
                     phases(i).vehicle.m0, phases(i).vehicle.gtow+1];                           % mass
phases(i).cbounds = [-20*pi/180, 20*pi/180; 0, 1; 0, 0]; % [min, max] for alpha, throttle, beta
phases(i).controls_no_match = [0, 1 ,0];

phases(i).tof = [120; 600];          % bounds for time of flight for the phase
phases(i).ode = @ode4;
phases(i).tstep= 1;                 % minimum timestep [seconds]


%% general constraints, cost function, settings and additional variables
param.c = {'c_acc', [6,6];
            'c_temp',3000;
            'c_maxq', 35e3};
param.objfun = 'max_payload';
% param.bias_obj_const = [1e0,1e0];

% first guess parameters---------------------------------------------------
param.fg = {'n_individuals',30;
            'distribution','random'};
%             'noise',1e-3;
%             'seed','case22'};

% additional extra optimisation variables for scaling of forces (L,D, T,
% mp) and second stage GTOW------------------------------------------------
% comment out the voices that are not needed, set equal bounds to multiply

param.Scale_S1_Engine = [7.6,8.2];
param.Scale_S2_Engine = [0.934,0.934];
% param.Scale_S1_Surf = [1,1];
% param.Scale_S2_Surf = [1,1];
% param.Scale_S2_GTOW = [1e3,1e3];

% parallel computing settings----------------------------------------------
param.parallel= 4;       % 0=force serial, 1=force parallel, other=manual
% param.debug = [];        % comment out for full optimisation and plotting
% no fly zone avoidance circles--------------------------------------------
% 1, 2 = lat lon of center of the circle, 3 = radius, 4 = max altitude where the limit apllies
% param.no_fly = [...
%                 43.1*pi/180,  17.1*pi/180, 2.30e6, 50e3;...  % europe
%                 65*pi/180,   -19.2*pi/180, 0.27e6, 50e3;...  % iceland
%                 64.2*pi/180,  27.2*pi/180, 0.80e6, 50e3;...  % norway
%                 18.5*pi/180,  -0.2*pi/180, 1.88e6, 50e3;...  % africa
%                 ]; 
            
%% load airport database for landing
% param.landing_zones = LoadAirports;
       
%% Run optimization routine -----------------------------------------------
run_TROPICO(param, phases)