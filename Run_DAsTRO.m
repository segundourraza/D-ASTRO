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

outputName ="test1";

plotOption = true;
saveOption = false;

%---------------------- AEROCAPTURE CORRDIOR METHOD ----------------------%
%   [1] = Single Re-entry trajectory analysis
%   [2] = Compute Aerocapture corridor SERIAL implementation
%   [3] = Compute Aerocapture corridor PARALLEL implementation
AerocaptureCorridor = 2;

gamma_range = -9.6494*pi/180;
BC_range = linspace(3, 60, 7);

densityMode = 1;    % Used when AerocaptureCorridor = 1 and robustCorridor = false
%---------------------------- OPTIMISATION MODE --------------------------%
%   [0] = NO
%   [1] = YES
optimisation = 1;

% x0 = [59.0000	-1.0472];
% x0 = [40.6000	-0.1765];
% x0 = [59.0000	-0.0010];
% x0 = [3.0010	-1e-3];
% x0 = [39.6000	-0.1894];
% x0 = [37.8000	-0.1833];

weights = [ 0;   % V_umbrella
            0;   % Fuel for corrective manoeuvress (Dv)
            1;   % Heat load (Q)
            0];     % Peak heatting rate (qdot_max)

targetOrbit = [ 4620.5;  % Semi-major axis (km)
                70;      % Inclination (deg)
                0.05];      % Eccentricity

%---------------------------- SETUP PHASE --------------------------%

run caseSetup.m

[simualtionOptions,Gammas ,OrbitalParams,Trajectories,DesignDrivers, ...
    BCArray, OptiResults] = simulatorUserInterface(fileData, simInputs);

OptiResults.GammaOpt =  OptiResults.GammaOpt *180/pi;
OptiResults