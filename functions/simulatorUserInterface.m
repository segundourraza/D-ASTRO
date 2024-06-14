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

    if any(simInputs.AerocaptureCorridor == [2,3])
        fitname = 'poly5';
        if simInputs.completeCorridor
            fitlb = fit(simInputs.BC_range', Gammas(:,5), fitname);
            fitub = fit(simInputs.BC_range', Gammas(:,4), fitname);
        else
            fitlb = fit(simInputs.BC_range', Gammas(:,1), fitname);
            fitub = fit(simInputs.BC_range', Gammas(:,2), fitname);
        end
        simInputs.fits = {fitlb, fitub};
    end
end


if simInputs.Opti.optimisation && simInputs.AerocaptureCorridor ~= 1
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

    % Run optimal Trajectory
    [~, trajResults] = Trajectory3DSim([OptiResults.BCOpt,OptiResults.GammaOpt], simInputs.densityMode, simInputs);
    t = trajResults(:,1);
    tend =  t(end);
    tq = linspace(0, tend, simInputs.nPoints);
    Vq = interp1(t, trajResults(:,2:end), tq, 'pchip');
    Trajectories{1} = [tq', Vq];
end
%% PLOTTING
if simInputs.plot_option
    if simInputs.AerocaptureCorridor == 1
        lineCycle = {'-', '--', ':', '-.'};
        markerCycle = {'none','o','^','.', 's'};
        nCycle = 5;
        ms = 6;
        
        fig = figure();
        sgtitle("Trajectory Profiles")
        set(fig, 'defaultlinelinewidth',2)
        set(fig, 'defaultaxesfontsize', 14)
        set(gca,'TickLabelInterpreter','latex')
        set(fig,'defaultTextInterpreter','latex'); %trying to set the default

        subplot(2,2,1);
        hold on
        h = zeros(1, length(Trajectories));
        for i = 1:length(Trajectories)
            if i <=nCycle-1
                style = [lineCycle(i),markerCycle(1)];
            else
                rem = mod(i,nCycle)+2;
                style = [lineCycle(floor(i/(nCycle-1))),markerCycle(rem)];
                if markerCycle(rem) == "."; ms = 24; end
            end
            h(1,i) = plot(Trajectories{i}(:,1), (Trajectories{i}(:,2) - simInputs.R)/1e3, 'DisplayName',sprintf("(BC, FPA) = (%5.1f, %2.2f deg)", simInputs.BC_range(i), simInputs.gamma_range(i)*180/pi), LineStyle=style(1), Marker=style(2), Color='k', MarkerSize=ms);
        end
        hold off
        ylim([0,inf])
        ylabel("$h$ (km)")
        xlabel( "$t$ (s)")
        
        lh=legend(h,'location','southeast');
        set(lh,'position',[.7 .25 .1 .1]);
        grid on ;grid minor;

        subplot(2,2,2)
        hold on
        for i = 1:length(Trajectories)
            yyaxis right
            plot(Trajectories{i}(:,1), Trajectories{i}(:,end-1)/1e4)
            ylabel("$\dot{q}$ (W/cm$^2$)" )

            yyaxis left
            plot(Trajectories{i}(:,1), Trajectories{i}(:,end)/1e6)
            ylabel("$Q$ (MJ/m$^2$)" )
        end
        hold off
        xlabel("$t$ (s)")
        grid on ;grid minor;

        subplot(2,2,3)
        hold on
        for i = 1:length(Trajectories)
            if i <=nCycle-1
                style = [lineCycle(i),markerCycle(1)];
            else
                rem = mod(i,nCycle)+2;
                style = [lineCycle(floor(i/(nCycle-1))),markerCycle(rem)];
                if markerCycle(rem) == "."; ms = 24; end
            end
            plot(Trajectories{i}(:,3)/1e3, (Trajectories{i}(:,2) - simInputs.R)/1e3,  LineStyle=style(1), Marker=style(2), Color='k', MarkerSize=ms)
        end
        hold off
        xlabel("$v$ (kms$^{-1}$)")
        ylabel("$h$ (km)")
        grid on ;grid minor;

    elseif any(simInputs.AerocaptureCorridor == [2,3])
        fig = figure();
        set(fig, 'defaultlinelinewidth',2)
        set(fig, 'defaultaxesfontsize', 20)
        set(fig, 'defaulttextinterpreter','tex')
        axis square
        run plotAeroCorridor.m
        hold on

        BC = linspace(BCArray(1), BCArray(end), 100);
        plot(BC, simInputs.fits{1}(BC)*180/pi, '-', 'color', "#7E2F8E")
        plot(BC, simInputs.fits{2}(BC)*180/pi, '-', 'color', "#7E2F8E")

        if simInputs.Opti.optimisation
            plot(Optimalvalues(2), Optimalvalues(3)*180/pi, '.k','MarkerSize',32, 'DisplayName', 'Design Point')
            xline(Optimalvalues(2), '--k' , 'LineWidth',1.5, HandleVisibility='off')
            yline(Optimalvalues(3)*180/pi, '--k' , 'LineWidth',1.5, HandleVisibility='off')
        end
        hold off


        fig = figure();
        sgtitle("Trajectory Profiles")
        set(fig, 'defaultlinelinewidth',2)
        set(fig, 'defaultaxesfontsize', 14)
        set(gca,'TickLabelInterpreter','latex')
        set(fig,'defaultTextInterpreter','latex'); %trying to set the default

        subplot(2,2,1)
        plot(Trajectories{1}(:,1), (Trajectories{1}(:,2) - simInputs.R)/1e3)

        ylim([0,inf])
        ylabel("$h$ (km)")
        xlabel( "$t$ (s)")
        grid on ;grid minor;
        
        subplot(2,2,2)
        yyaxis right
        plot(Trajectories{1}(:,1), Trajectories{1}(:,end-1)/1e4)
        ylabel("$\dot{q}$ (W/cm$^2$)" )

        yyaxis left
        plot(Trajectories{1}(:,1), Trajectories{1}(:,end)/1e6)
        ylabel("$Q$ (MJ/m$^2$)" )
        hold off
        grid on ;grid minor;
        xlabel("$t$ (s)")

        subplot(2,2,3)
        plot(Trajectories{1}(:,3)/1e3, (Trajectories{1}(:,2) - simInputs.R)/1e3)
        ylim([0,inf])
        xlabel("$v$ (kms$^{-1}$)")
        ylabel("$h$ (km)")
        grid on ;grid minor;

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
