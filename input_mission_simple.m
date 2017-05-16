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

%% Rocket ascent
i = i+1;
phases(i).vehicle = load_vehicle('test');
phases(i).ne   = 2;             % number of elements per phase (for multi-shooting)
phases(i).nc   = 5;            % number of control nodes for each element within the phase
phases(i).x0 =          [15e3,  0.7*300, 0, pi/2, -5*pi/180, 55*pi/180, 100e3];
phases(i).x0_no_opt =   [1, 1, 1, 0, 0, 0, 1];                              % id of fixed starting values
phases(i).xf =          [25e3, 1e3, 0, pi/2, 0, 0, 50e3];
phases(i).xf_fixed =    [1, 0, 0, 0, 0, 0, 0];                              % id of target  value
phases(i).xbounds = [1e3, 100e3;...            % altitude
                     100, 8e3; ...           % velocity
                     -pi/2, pi/2; 0, 2*pi; ...                                  % FPA, heading
                     -pi/2, pi/2; -pi, pi; ...       % lat, lon
                     1, 100e3];                           % mass
phases(i).cbounds = [-10*pi/180, 30*pi/180; 0, 1; 0, 0]; % [min, max] for alpha, throttle, beta
phases(i).fg = [(10*pi/180), 1, 0];     % starting guess for controls, can be any number of points > 2 for each control variable
% phases(i).controls_no_match = [0, 1 ,0];
phases(i).tof = [10; 120];          % bounds for time of flight for the phase

phases(i).ode = @ode4;
phases(i).odet = phases(i).nc*2;
phases(i).tstep= 1;                 % minimum timestep [seconds]
phases(i).atmo = @atmo_ext_ISA;
phases(i).dynamics = @dyn_ECEF2;
phases(i).prop = @prop_rocket_450s_1MN;
phases(i).aero = @aero_simplerocket;
phases(i).thermal = @temp_convradbalance;
phases(i).inte = @interp_pch_mex;
phases(i).dist = @linspace; %@dist_chebyshev_mod;

% phases(i).ceq = {'ceq_circ',0};

%% general constraints, cost function, settings and additional variables
param.c = {'c_acc', [6,6];
            'c_temp',3000;
            'c_maxq', 20e3};
param.objfun = 'min_fuel';
param.bias_obj_const = [1e0,1e0];
% param.objfun = 'max_c3';

% first guess parameters---------------------------------------------------
% field 1 is approach, field 2 is timestep multiplier
% param.fg = {'lhs_feasibility_fmincon',1};   
param.fg = {'lhs_fmincon',1};       
% param.fg = {'none'};
% param.fg = {'ga',3};    
% param.fg = {'load','success_20170503T122046.mat'};

% only for demonstration purposes:
% param.fg = {'load','success_20170120T164142.mat'};  % test result - ascent
% param.fg = {'load','fail_20170123T050714.mat'};  % test result - descent

% additional extra optimisation variables for scaling of forces (L,D, T,
% mp) and second stage GTOW------------------------------------------------
% comment out the voices that are not needed, set equal bounds to multiply

param.Scale_S1_Engine = [3,3];
% param.Scale_S2_Engine = [1,1];
% param.Scale_S1_Surf = [1,1];
% param.Scale_S2_Surf = [1,1];
% param.Scale_S2_GTOW = [1e3,1e3];

% parallel computing settings----------------------------------------------
param.parallel= 1;       % 0=force serial, 1=force parallel, other=manual
% param.debug = [];        % comment out for full optimisation and plotting
% no fly zone avoidance circles--------------------------------------------
% 1, 2 = lat lon of center of the circle, 3 = radius, 4 = max altitude where the limit apllies
param.no_fly = [...
%                 43.1*pi/180,  17.1*pi/180, 2.30e6, 50e3;...  % europe
%                 65*pi/180,   -19.2*pi/180, 0.27e6, 50e3;...  % iceland
%                 64.2*pi/180,  27.2*pi/180, 0.80e6, 50e3;...  % norway
%                 18.5*pi/180,  -0.2*pi/180, 1.88e6, 50e3;...  % africa
                ]; 
            
%% load airport database for landing
% param.landing_zones = LoadAirports;
       
%% Run optimization routine -----------------------------------------------
run_TROPICO(param, phases)