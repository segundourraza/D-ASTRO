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

if exist('outputName','var') == 1
    fileData.outputName = outputName;
else
    fileData.outputName = replace(string(datetime('now','Format','hh:mm:ss')),':','-');
end

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


if exist('planet','var') == 0
    planet = "mars";
end

if exist('defaultMarsModel','var') == 0
    defaultMarsModel = 1;
end


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

if defaultMarsModel && lower(planet) == "mars"
    addpath("atmospheric models\Mars\")
    simInputs = inputSim();
    run 'atmospheric models'\Mars\atmoMarsTables.m  
    simInputs.atmoModel = @(h, denMod, full_analysis) atmoModelMars(h, denMod, full_analysis, simInputs.models);
elseif lower(planet) == "earth"
    simInputs = inputSim_Earth();
    simInputs.atmoModel = @(h, denMod, full_analysis) atmoModelEarth(h, denMod);
end


simInputs.plot_option = plotOption;
simInputs.save_option = saveOption;
%% Aerocapture corridor paramteres
simInputs.AerocaptureCorridor = AerocaptureCorridor;
if simInputs.AerocaptureCorridor ==1
    simInputs.gamma_range = gamma_range;
end
simInputs.BC_range =  BC_range;
simInputs.densityMode = densityMode;

%% OPTIMIZATION PARAMETERS
if optimisation
    % Initial condition
    if exist('x0','var')
        simInputs.Opti.x0 = x0;
    end
    % Optimisation paramteres
    simInputs.Opti.optimisation = optimisation;
    simInputs.Opti.weights = [weights', 0, 0, 0];
    simInputs.Opti.target = [0, 0, 0, 0, 0, 0,0];
    
    % Operational target orbit
    simInputs.Opti.BC_range = BC_range;
    simInputs.Opti.a_orb_target = targetOrbit(1)*1e3;
    simInputs.Opti.i_orb_target = targetOrbit(2)*pi/180;
    simInputs.Opti.e_orb_target = targetOrbit(3);
end
