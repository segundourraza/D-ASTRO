%-------------------------------------------------------------------------%
%                    Determination of Aerocapture Successful              %
%                         manouvre (D-ASTRO) and Robust                   %
%                            Optimisation algorithm                       %
%                                                                         %
%                              Segundo Urraza Atue                        %
%                         Imperial College London, 2023                   %
%-------------------------------------------------------------------------%

run preamble.m

simInputs = inputSim();

%---------------------- AEROCAPTURE CORRDIOR METHOD ----------------------%
%   [0] = Optimise exisitng corridor
%   [1] = Single Re-entry trajectory analysis
%   [2] = Compute Aerocapture corridor SERIAL implementation
%   [3] = Compute Aerocapturetest5 corridor PARALLEL implementation
%
simInputs.AerocaptureCorridor = 3;
simInputs.AerocapFiles = "savedExperiment_step01V3";

%---------------------------- OPTIMISATION MODE --------------------------%
%   [0] = NO 
%   [1] = YES

simInputs.Opti.optimisation = 1;


gamma_range = -9.6494*pi/180;
BC_range = 5:1:simInputs.BCmax;
BC_range = 3:10:60;
BC_range = linspace(3, 60, 7);

% BC_range = linspace(BC_range(1), BC_range(end), 100);
% BC_range = 11.589;

simInputs.gamma_range = gamma_range; 
simInputs.BC_range =  BC_range;
simInputs.Opti.BC_range = BC_range;
simInputs.Opti.a_orb_target = 4.6205e+06;
simInputs.Opti.i_orb_target = 70*pi/180;
simInputs.Opti.e_orb_target = 0.05;

% x = [V_umbrella, DV_total, Q, q_dot_max, Tw, P,Dgamma];
simInputs.Opti.weights = [1, 1, 1, 1, 0, 0, 0];
simInputs.Opti.Target = [0, 0, 0, 0, 0, 0, 0];

% MUST PROVIDE ATMOPSHERIC A FUNCTION TO MODEL THE ATMOPSHERE OF THE
% TARGET PLANET. 
% FUNCTION MUST HAVE NAME "marsatm" AND  MUST HAVE THE FOLLOWING FORMAT:
%
% [T, P, rho] = marsatm(h, option, full_analysis)
%   INPUTS:     -h: Altitest3tude (m)
%               -option: Determines whether the average density is used or
%               an upper/lower bound based from Mars GRAM deviations.
%                       [1] = mean atmopshere (default)
%                       [2] = lower limit of atmopshere
%                       [3] = higher limit of atmopshere
%               -full_analysis: Option for which analysis is to be done.
%   OUTPUTS:    -T: Temperature (K)
%               -P: pressure (Pa)
%               -rho: density (Kg/m3)
run loadAtmModel.m
simInputs.models = models;

[simualtionOptions,Gammas ,OrbitalParams,Trajectories,DesignDrivers, ...
    BCArray, Optimalvalues, OptiResults] = simulatorUserInterface(fileData, simInputs);

OptiResults.GammaOpt =  OptiResults.GammaOpt *180/pi;
OptiResults