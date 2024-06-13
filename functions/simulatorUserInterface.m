function [simualtionOptions,Gammas ,OrbitalParams,Trajectories,DesignDrivers, BCArray, Optimalvalues, OptiResults] = simulatorUserInterface(fileData, simInputs)
simInputs.Opti.weights = (simInputs.Opti.weights - min(simInputs.Opti.weights))./(max(simInputs.Opti.weights) - min(simInputs.Opti.weights));

simualtionOptions = [];
Gammas = [];
OrbitalParams = [];
Trajectories = [];
DesignDrivers = [];
BCArray =  [];
Optimalvalues = [];
OptiResults = [];

fileData.outputName = input('Save simualtion results as: ', 's');
% 
% startTime = now();
% startTimeS = datestr(now,'HH:MM:SS.FFF');
% disp("Start time: "+startTimeS);
% 
% if simInputs.AerocaptureCorridor ~=0
%     if simInputs.AerocaptureCorridor ~= 1
%         fprintf("\nCOMPUTING AEROCAPTURE CORRIDOR...\n")
%     else
%         simInputs.density_analysis = 0;
%     end
% %     [Gammas, OrbitalParams, Trajectories, DesignDrivers] = AtmTrajSimV4(simInputs);
%     [Gammas, OrbitalParams, Trajectories, DesignDrivers] = AtmTrajSimV5(simInputs);
%     var_names = {'BC','Flight path angle','Radius body', 'Delta V[Km/s]','Heat Load[MJ/cm2]', 'Max Heat transfer[W/cm2]', 'Max Wall Temp.[K]', 'Delta Gamma [deg]' };
%     T1 = table(simInputs.BC_range,simInputs.gamma_range*180/pi, DesignDrivers(1), DesignDrivers(2), DesignDrivers(3)/1e6*1e4, DesignDrivers(4)/1e4, DesignDrivers(5), DesignDrivers(6)*180/pi, 'VariableNames', var_names);
%     disp(T1)
%     BCArray = simInputs.BC_range;
%     fprintf("\n")
% end
% simualtionOptions = simInputs;
% 
% if simInputs.AerocaptureCorridor==1 && length(simInputs.BC_range) == 1
% 
%     fig = figure();
%     set(fig, 'defaultlinelinewidth',2)
%     set(fig, 'defaultaxesfontsize', 14)
%     set(fig, 'defaulttextinterpreter','tex')
%     t = Trajectories{1}(:,1);
%     h = Trajectories{1}(:,2) - simInputs.R;
%     v = Trajectories{1}(:,3);
%     plot(t, h/1e3)
% 
% 
%     time_atm_pass = Trajectories{1}(end,1); % ATmopsheric pass time
%     distance_atm_pass = trapz(Trajectories{1}(:,1),Trajectories{1}(:,3)); % Atmospheric pass distance% Thermal model
%     var_names = {'Flight path angle [deg]','Semi-major axis[Km]', 'eccentricity','Inclination' 'True anomaly (deg)', 'Atmospheric time pass (s)', 'Distance atmospheric pass (Km)'};
%     T1 = table(simInputs.gamma_range*180/pi,OrbitalParams(1)/1e3,OrbitalParams(2),OrbitalParams(3)*180/pi,OrbitalParams(4)*180/pi,time_atm_pass,distance_atm_pass/1e3,'VariableNames', var_names);
%     %     T1 = varfun(@(x) num2str(x, ['%' sprintf('.%de', 2)]), T1);
%     %     T1.Properties.VariableNames = var_names;
%     disp(T1)
%     % Save Results
%     if ~isfolder(fileData.OutputPath)
%         mkdir(fileData.OutputPath);
%     end
% 
%     saveFilePath = fullfile(pwd,fileData.OutputPath,"savedExperiment_" +fileData.outputName);
%     save(saveFilePath,'simualtionOptions','Gammas','OrbitalParams','Trajectories', 'DesignDrivers', 'BCArray');
%     disp("Results saved in path: " + saveFilePath)
%     fprintf("\n")
% elseif any(simInputs.AerocaptureCorridor == [2,3]) && ~simInputs.Opti.optimisation
%     fig = figure();
%     set(fig, 'defaultlinelinewidth',2)
%     set(fig, 'defaultaxesfontsize', 20)
%     set(fig, 'defaulttextinterpreter','tex')
%     axis square
%     run AerocaptureCorridor.m
% 
%     % Save Results
%     if ~isfolder(fileData.OutputPath)
%         mkdir(fileData.OutputPath);
%     end
% 
%     saveFilePath = fullfile(pwd,fileData.OutputPath,"savedExperiment_" +fileData.outputName);
%     save(saveFilePath,'simualtionOptions','Gammas','OrbitalParams','Trajectories', 'DesignDrivers', 'BCArray');
%     disp("Results saved in path: " + saveFilePath)
%     fprintf("\n")
% elseif any(simInputs.AerocaptureCorridor == [2,3]) && simInputs.Opti.optimisation
%     fprintf("\nSTART OPTIMISATION....\n")
% %     [Optimalvalues, OptiResults] = optimise(BCArray, Gammas, simInputs);
%     [Optimalvalues, OptiResults] = optimiseGD(BCArray, Gammas, DesignDrivers, simInputs);
% 
%     fig = figure();
%     set(fig, 'defaultlinelinewidth',2)
%     set(fig, 'defaultaxesfontsize', 14)
%     set(fig, 'defaulttextinterpreter','tex')
% 
%     run AerocaptureCorridor.m
% 
%     figure(fig)
%     hold on
%     plot(Optimalvalues(2), Optimalvalues(3)*180/pi, '.k','MarkerSize',32, 'DisplayName', 'Design Point')
%     xline(Optimalvalues(2), '--k' , 'LineWidth',1.5)
%     yline(Optimalvalues(3)*180/pi, '--k' , 'LineWidth',1.5)
%     hold off
% 
%     legend off
% 
% 
%     % Save Results
%     if ~isfolder(fileData.OutputPath)
%         mkdir(fileData.OutputPath);
%     end
%     saveFilePath = fullfile(pwd,fileData.OutputPath,"savedExperiment_" +fileData.outputName);
%     save(saveFilePath,'simualtionOptions','Gammas','OrbitalParams','Trajectories','DesignDrivers','BCArray', "Optimalvalues", "OptiResults");
%     disp("Results saved in path: " + saveFilePath)
%     fprintf("\n")
% else
%     load("Output files\" + simInputs.AerocapFiles);
% 
%     Opti = simInputs.Opti;
%     simInputs = simualtionOptions;
%     simInputs.verbose = 1;
%     simInputs.Opti = Opti;
%     fprintf("\nSTART OPTIMISATION....\n")
% %     [Optimalvalues, OptiResults] = optimise(BCArray, Gammas, simInputs);
%     tic    
%     for i = 1:1
%         [Optimalvalues, OptiResults] = optimiseGD(BCArray, Gammas, DesignDrivers, simInputs);
%     end
%     toc
%     fig = figure();
%     set(fig, 'defaultlinelinewidth',2)
%     set(fig, 'defaultaxesfontsize', 14)
%     set(fig, 'defaulttextinterpreter','tex')
% 
%     run AerocaptureCorridor.m
% 
%     figure(fig)
%     hold on
%     plot(Optimalvalues(2), Optimalvalues(3)*180/pi, '.k','MarkerSize',32, 'DisplayName', 'Design Point')
%     xline(Optimalvalues(2), '--k' , 'LineWidth',1.5)
%     yline(Optimalvalues(3)*180/pi, '--k' , 'LineWidth',1.5)
%     hold off
% 
%     legend off
% 
%     % Save Results
%     if ~isfolder(fileData.OutputPath)
%         mkdir(fileData.OutputPath);
%     end
%     saveFilePath = fullfile(pwd,fileData.OutputPath,"savedExperiment_" +fileData.outputName);
%     save(saveFilePath,'simualtionOptions','Gammas','OrbitalParams','Trajectories','DesignDrivers','BCArray', "Optimalvalues", "OptiResults");
%     disp("Results saved in path: " + saveFilePath)
%     fprintf("\n")
% end
% 
% elapsedTimeSeconds = (now() - startTime)*1e5;
% endTimes = datestr(now,'HH:MM:SS.FFF');
% fprintf('\n')
% disp("End time: "+endTimes);
% disp("Elapsed time: " + elapsedTimeSeconds);
% 
% clearvars -except simualtionOptions Gammas OrbitalParams Trajectories DesignDrivers BCArray Optimalvalues OptiResults
