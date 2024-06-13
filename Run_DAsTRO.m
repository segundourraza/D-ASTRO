%-------------------------------------------------------------------------%
%                    Determination of Aerocapture Successful              %
%                         manouvre (D-ASTRO) and Robust                   %
%                            Optimisation algorithm                       %
%                                                                         %
%                              Segundo Urraza Atue                        %
%                         Imperial College London, 2023                   %
%-------------------------------------------------------------------------%
close all
clear
clc

outputName ="";

plotOption = false;
saveOption = false;

%---------------------- AEROCAPTURE CORRDIOR METHOD ----------------------%
%   [0] = Optimise exisitng corridor
%   [1] = Single Re-entry trajectory analysis
%   [2] = Compute Aerocapture corridor SERIAL implementation
%   [3] = Compute Aerocapturetest5 corridor PARALLEL implementation
%
AerocaptureCorridor = 1;
AerocapFiles = "savedExperiment_step01V3";

gamma_range = -9.6494*pi/180;
BC_range = linspace(3, 60, 7);

% Atmopsheric model options
% Use dafult Martian model used in D-ASTRO publication?
planet = 'mars';
defaultMarsModel = true;
defaultEarthModel = true;

densityMode = 3;    % Only used when AerocaptureCorridor = 1
%---------------------------- OPTIMISATION MODE --------------------------%
%   [0] = NO 
%   [1] = YES

optimisation = 1;

weights = [ 1,...   % V_umbrella  
            1,...   % Fuel for corrective manoeuvress (Dv)
            1,...   % Heat load (Q)
            1];     % Peak heatting rate (qdot_max)

targetOrbit = [ 4620.5,...  % Semi-major axis (km)
                70,...      % Inclination (deg)
                0.05];      % Eccentricity

%---------------------------- SETUP PHASE --------------------------%

simInputs = inputSim();
run caseSetup.m

[simualtionOptions,Gammas ,OrbitalParams,Trajectories,DesignDrivers, ...
    BCArray, Optimalvalues, OptiResults] = simulatorUserInterface(fileData, simInputs);
% 
% OptiResults.GammaOpt =  OptiResults.GammaOpt *180/pi;
% OptiResults