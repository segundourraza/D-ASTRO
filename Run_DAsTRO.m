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

plotOption = false;
saveOption = false;

%---------------------- AEROCAPTURE CORRDIOR METHOD ----------------------%
%   [0] = Optimise exisitng corridor
%   [1] = Single Re-entry trajectory analysis
%   [2] = Compute Aerocapture corridor SERIAL implementation
%   [3] = Compute Aerocapture corridor PARALLEL implementation
AerocaptureCorridor = 2;

AerocapFiles = "savedExperiment_step01V3";

gamma_range = -9.6494*pi/180;
BC_range = linspace(3, 60, 7);

if AerocaptureCorridor == 1
    xopt = [13.03 -9.756;
            48.46 -10.487;
            3.00 -8.7822;
            3.0 -8.780;
            3.000 -8.780];
    BC_range = xopt(:,1);
    gamma_range = xopt(:,2)*pi/180;


    gamma_range = linspace(-11.5, -10, 7)*pi/180;
    BC_range = ones([1,length(gamma_range)])*20;
end

densityMode = 1;    % Used when AerocaptureCorridor = 1 and robustCorridor = false
%---------------------------- OPTIMISATION MODE --------------------------%
%   [0] = NO
%   [1] = YES
optimisation = 1;

weights = [ 1;   % V_umbrella
            1;   % Fuel for corrective manoeuvress (Dv)
            1;   % Heat load (Q)
            1];     % Peak heatting rate (qdot_max)

targetOrbit = [ 4620.5;  % Semi-major axis (km)
                70;      % Inclination (deg)
                0.05];      % Eccentricity

%---------------------------- SETUP PHASE --------------------------%

simInputs = inputSim();
run caseSetup.m

[simualtionOptions,Gammas ,OrbitalParams,Trajectories,DesignDrivers, ...
    BCArray, Optimalvalues, OptiResults] = simulatorUserInterface(fileData, simInputs);

OptiResults.GammaOpt =  OptiResults.GammaOpt *180/pi;
OptiResults