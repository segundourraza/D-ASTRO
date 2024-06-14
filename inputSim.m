function simInputs = inputSim()

%--------------------------- PLANETARY PARAMETERS ------------------------%

simInputs.R = 3389.5e3;                     % Radius of target planet (m)
simInputs.mu = 4.282837e13;                 % Gravitational parameter of target planet (m3 s−2)
simInputs.Omega = 7.0830e-05;               % Angular frequency of target planet (rad s-1)

simInputs.k = 1.83e-4;                      % Sutton-Graves constant
simInputs.sigma = 5.67e-8;                  % Boltzmann constant
simInputs.C = 2.35e4;                       % Sutton-Tauber constant
simInputs.f_V =@(V) tauberSuttonFunction(V);    % empirical function of velcoity for radiative heating

%-------------------------- SPACECRAFT PARAMETERS ------------------------%

simInputs.SC.m = 400;           % Mass of spacecraft (kg)
simInputs.SC.L = 2.6;           % Diameter of fairing (m)
simInputs.SC.Tw = 290;          % Starting wall tempertaure of aeroshell (k)
simInputs.SC.CD = 1.6;          % Drag coefficient of aeroshell
simInputs.SC.L_D = 0.2;         % Lift-to-DRag ratio of spacecraft

simInputs.SC.phi = deg2rad(70);
simInputs.SC.RN_RB = 0.5;       % Ratio of nose radius to body radius;

simInputs.SC.maxFairingRadius = 1.3;
simInputs.BCmax = simInputs.SC.m/(simInputs.SC.CD*pi*simInputs.SC.maxFairingRadius^2);

simInputs.SC.epsilon = 0.9;             % Emissivity of TPS

%--------------------- INSERTION TRAJECTORY PARAMETERS -------------------%
simInputs.Dgamma = deg2rad(0.2);

simInputs.Vinf = 3.5;                  % Hyperbolit excess velocity (km s-1)
% simInputs.Vinf = 6.0;                  % Hyperbolit excess velocity (km s-1)
simInputs.h_AI = 125e3;                % Altitude of atmopsheric interface (AI)

simInputs.R0 = simInputs.h_AI + simInputs.R;
simInputs.V0 = sqrt(2*simInputs.mu/simInputs.R0 + (simInputs.Vinf*1e3)^2);
simInputs.delta0 = deg2rad(0.5798);    % Initial Longitude, +ve to east (rad)
simInputs.lambda0 = deg2rad(34.49);    % Initial Latitude, +ve northern hemisphere (rad)
simInputs.chi0 = deg2rad(90-108.24);   % Heading angle (rad)

%---------------------------- GENERAL OPTIONS ----------------------------%

simInputs.robustCorridor = false;
simInputs.nPoints = 100;

simInputs.Opt    = odeset('Events', @(t,y) myEventLimit(t,y,[simInputs.R, simInputs.h_AI]),...
                            'RelTol',1e-6,'AbsTol',1e-6);
simInputs.tspan = [0, 1e3];

simInputs.thermalModel = 1;     % [1]:  Sutton-graves

simInputs.VolumeMethod = 3;     % [1]:  Stowed deployable aroshell volume
                                % [2]:  Deployed aeroshell volume
                                % [3]:  Deployed aeroshell body radius (RECOMMENDED)
%------------------- BINARY SEARCH ALGORITHMS OPTIONS --------------------%

simInputs.gamma_min = -90*pi/180;
simInputs.gamma_max = -0*pi/180;

simInputs.tol = 1e-8;       % FPA binary search algorithm tolerence
simInputs.MaxIterations = 1e4;

%-------------------- OPTIMISATION ALGORITHM OPTIONS ---------------------%

DispOpt ='iter';
simInputs.Opti.options = optimoptions('fmincon',...
                                      'Display',DispOpt,...
                                      'MaxIterations',8000, ...
                                      'PlotFcn','optimplotfvalconstr',...
                                      'FunctionTolerance',1e-16, 'OptimalityTolerance',1e-16);
                        
simInputs.Opti.costFunc = 1;
simInputs.Opti.delta = 0.1;
simInputs.verbose = 0;

end