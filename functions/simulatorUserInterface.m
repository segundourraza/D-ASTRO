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


startTime = datetime('now','Format','hh:mm:ss.SSS');
fprintf("Start time: \t %s\n", string(startTime));

if fileData.outputName == ""
    fileData.outputName = replace(string(datetime('now','Format','hh:mm:ss')),':','-');
end

if simInputs.AerocaptureCorridor ~= 0
    if simInputs.AerocaptureCorridor ~= 1
        fprintf("\nCOMPUTING AEROCAPTURE CORRIDOR...\n")
    else
        if length(simInputs.BC_range) ~= length(simInputs.gamma_range)
            if length(simInputs.gamma_range) == 1
                simInputs.gamma_range = ones(length(simInputs.BC_range))*simInputs.gamma_range;
            else
                error("ERROR: List of gammas must be equal to length of BC list or of length (1). Currently of length (" + length(simInputs.gamma_range) + ").")
            end
        end
        simInputs.density_analysis = 0;
    end
    [Gammas, OrbitalParams, Trajectories, DesignDrivers] = aeroCorridor(simInputs);
    BCArray = simInputs.BC_range;
    fprintf("\n")
end


simualtionOptions = simInputs;
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


% Save Results
if simInputs.save_option
    if ~isfolder(fileData.OutputPath)
        mkdir(fileData.OutputPath);
    end
    
    saveFilePath = fullfile(pwd,fileData.OutputPath,"savedExperiment_" +fileData.outputName);
    
    if simInputs.AerocaptureCorridor == 1
        for i = 1:length(Trajectories)
            T = array2table(Trajectories{i});
            T.Properties.VariableNames(1:9) = {'t(s)', 'r(m)', 'v (m/s)', 'FPA (rad)', 'Longitude (rad)', 'Latitude (rad)', 'Heading angle (rad)', 'Heat transfer (W/m2)', 'Heat load (J/m2)'};
            writetable(T, saveFilePath+ "_Trajectory" + i +".csv")
            
        end
    end    
    
    save(saveFilePath,'simualtionOptions','Gammas','OrbitalParams','Trajectories','DesignDrivers','BCArray', "Optimalvalues", "OptiResults");
    disp("Results saved in path: " + saveFilePath)
    fprintf("\n")
end
endTime = datetime('now','Format','hh:mm:ss.SSS');
elapsedTimeSeconds = endTime - startTime;
elapsedTimeSeconds.Format = 'hh:mm:ss.SSS';
fprintf('\n')
fprintf("End time: \t\t %s\n", string(endTime));
fprintf("Elapsed time: \t %s\n", string(elapsedTimeSeconds));

% clearvars -except simualtionOptions Gammas OrbitalParams Trajectories DesignDrivers BCArray Optimalvalues OptiResults
