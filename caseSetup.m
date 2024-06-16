%% HOUSEKEEPING
restoredefaultpath;

if any(computer== ["PCWIN64","MACI64"])  % Windows or Mac OS
    fileData.computerName = getenv('COMPUTERNAME');
    % Use defaults
    fileData.OutputPath = fullfile("Output files"); % Output directory
    fileData.WorkingPath = fullfile(pwd + "/"); % Working directory
elseif computer =="GLNXA64" % Linux OS
    fileData.computerName = getenv('COMPUTERNAME');
    %Use defaults
    fileData.OutputPath = fullfile("Output file"); % Output directory
    fileData.WorkingPath = fullfile(pwd + "\"); % Working directory
end
fileData.outputName = outputName;

% Check for dependencies
v = ver;
[installedToolboxes{1:length(v)}] = deal(v.Name);
if AerocaptureCorridor == 3 && ~ismember('Parallel Computing Toolbox', installedToolboxes)
    warning("'Parallel Computing Toolbox' is not installed. Switching to SERIAL implementation.")
    AerocaptureCorridor =2;
end

if optimisation && ~ismember('Optimization Toolbox', installedToolboxes)
    error("'Optimization Toolbox' is not installed. Only trajectory analysis tool can be executed (AerocaptureCorridor = 1).")
end

% Add required paths
addpath("atmospheric models\")
addpath("functions\")
addpath("functions\trajectory\")
addpath("functions\optimization\")


simInputs.plot_option = plotOption;
simInputs.save_option = saveOption;
%% Aerocapture corridor paramteres
simInputs.AerocaptureCorridor = AerocaptureCorridor;
simInputs.AerocapFiles = AerocapFiles;
simInputs.gamma_range = gamma_range;
simInputs.BC_range =  BC_range;
simInputs.densityMode = densityMode;

%% OPTIMIZATION PARAMETERS

% Optimisation paramteres
simInputs.Opti.optimisation = optimisation;
simInputs.Opti.weights = [weights', 0, 0, 0];
simInputs.Opti.target = [0, 0, 0, 0, 0, 0,0];

% Operational target orbit
simInputs.Opti.BC_range = BC_range;
simInputs.Opti.a_orb_target = targetOrbit(1)*1e3;
simInputs.Opti.i_orb_target = targetOrbit(2)*pi/180;
simInputs.Opti.e_orb_target = targetOrbit(3);

%% ATMOPSHERIC MODELS

% MUST PROVIDE ATMOPSHERIC FUNCTION TO MODEL THE ATMOPSHERE OF THE
% TARGET PLANET. FUNCTION MUST HAVE NAME "atmoModel" AND  
% MUST HAVE THE FOLLOWING FORMAT:

% [T, P, rho] = atmoModel(h, option, full_analysis)
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

if simInputs.defaultMarsModel && lower(simInputs.planet) == "mars"
    addpath("atmospheric models\Mars\")
    run 'atmospheric models'\Mars\atmoMarsTables.m  
    simInputs.atmoModel = @(h, denMod, full_analysis) atmoModelMars(h, denMod, full_analysis, simInputs.models);
elseif lower(simInputs.planet) == "earth"
    simInputs.atmoModel = @(h, denMod, full_analysis) atmoModelEarth(h, denMod);
    % error('ERROR: Earth atmosphere not yet coded!')
end