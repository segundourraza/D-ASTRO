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
    if simInputs.AerocaptureCorridor == 1
        if length(simInputs.BC_range) ~= length(simInputs.gamma_range)
            if length(simInputs.gamma_range) == 1
                simInputs.gamma_range = ones(length(simInputs.BC_range))*simInputs.gamma_range;
            else
                error("ERROR: List of gammas must be equal to length of BC list or of length (1). Currently of length (" + length(simInputs.gamma_range) + ").")
            end
        end
        simInputs.density_analysis = 0;
    else
        fprintf("\nCOMPUTING AEROCAPTURE CORRIDOR...\n")
    end
    [Gammas, OrbitalParams, Trajectories, DesignDrivers] = aeroCorridor(simInputs);
    BCArray = simInputs.BC_range;
end


if simInputs.Opti.optimisation
    fprintf("\nSTART OPTIMISATION....\n")
    if simInputs.completeCorridor
        nTheta = 7;
        temp = zeros([size(DesignDrivers,1),nTheta*2]);
        temp(:,1:nTheta) = DesignDrivers(:,nTheta*(5-1)+1:nTheta*5);
        temp(:,nTheta+1:nTheta*2) = DesignDrivers(:,nTheta*(4-1)+1:nTheta*4);
        DesignDrivers = temp;
        [Optimalvalues, OptiResults] = optimiseGD(BCArray, Gammas(:,[5,4]), DesignDrivers, simInputs);
    else
        [Optimalvalues, OptiResults] = optimiseGD(BCArray, Gammas, DesignDrivers, simInputs);
    end
end
%% PLOTTING
if simInputs.plot_option
    if simInputs.AerocaptureCorridor == 1

        lineCycle = {'-', '--', ':', '-.'};
        markerCycle = {'none','o','^','.', 's'};
        nCycle = 5;
        ms = 6;
        fig = figure();
        set(fig, 'defaultlinelinewidth',2)
        set(fig, 'defaultaxesfontsize', 14)
        set(fig, 'defaulttextinterpreter','tex')

        subplot(1,2,1)
        hold on
        for i = 1:length(Trajectories)
            if i <=nCycle-1
                style = [lineCycle(i),markerCycle(1)];
            else
                rem = mod(i,nCycle)+2;
                style = [lineCycle(floor(i/(nCycle-1))),markerCycle(rem)];
                if markerCycle(rem) == "."; ms = 24; end
            end
            plot(Trajectories{i}(:,1), (Trajectories{i}(:,2) - simInputs.R)/1e3, 'DisplayName',sprintf("BC = %6.2f", simInputs.BC_range(i)), LineStyle=style(1), Marker=style(2), Color='k', MarkerSize=ms)
        end
        hold off
        ylim([0,inf])
        ylabel("$h (km)$", Interpreter="latex")
        xlabel( "$t (s)$", Interpreter="latex")
        legend('Location','best')
        subplot(1,2,2)
        hold on
        for i = 1:length(Trajectories)
            yyaxis right
            plot(Trajectories{i}(:,1), Trajectories{i}(:,end-1)/1e4,'DisplayName',sprintf("BC = %6.2f", simInputs.BC_range(i)))
            ylabel("$\dot{q}_{max} (W/cm^2)$", Interpreter="latex" )

            yyaxis left
            plot(Trajectories{i}(:,1), Trajectories{i}(:,end)/1e6)
            ylabel("$Q (MJ/m^2)$", Interpreter="latex" )
        end
        hold off
        xlabel("$t (s)$", Interpreter="latex")
        legend('Location','best')

        fig = figure();
        set(fig, 'defaultlinelinewidth', 2);
        hold on
        for i = 1:length(Trajectories)
            yyaxis right
            plot(Trajectories{i}(:,1), Trajectories{i}(:,end-1)/1e4)
            ylabel("$\dot{q}_{max} (W/cm^2)$", Interpreter="latex" )

            yyaxis left
            plot(Trajectories{i}(:,1), Trajectories{i}(:,end)/1e6)
            ylabel("$Q (MJ/m^2)$", Interpreter="latex" )
        end
        hold off
        xlabel("$t (s)$", Interpreter="latex")

    elseif any(simInputs.AerocaptureCorridor == [2,3])
        fig = figure();
        set(fig, 'defaultlinelinewidth',2)
        set(fig, 'defaultaxesfontsize', 20)
        set(fig, 'defaulttextinterpreter','tex')
        axis square
        run plotAeroCorridor.m

        hold on
        if simInputs.Opti.optimisation
            plot(Optimalvalues(2), Optimalvalues(3)*180/pi, '.k','MarkerSize',32, 'DisplayName', 'Design Point')
            xline(Optimalvalues(2), '--k' , 'LineWidth',1.5, HandleVisibility='off')
            yline(Optimalvalues(3)*180/pi, '--k' , 'LineWidth',1.5, HandleVisibility='off')
        end
        hold off

    end
end

%% Save Results
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
