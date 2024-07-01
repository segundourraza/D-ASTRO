clc
clear
close all
cd ("..\")

plotOption = false;
saveOption = false;

AerocaptureCorridor = 2;
densityMode = 3;

simInputs = inputSim_Earth();
BC_range = linspace(3, 60, 7);

optimisation = 1;
weights = [1, 1, 1, 1]';
targetOrbit = [6371 + 400, 70, 0.05];
planet = 'earth';
run caseSetup.m
cd("Performance\")


addpath("Multi-objective optimization\")
addpath("optimizers no fits\")

tic
[Gammas, OrbitalParams, Trajectories, DesignDrivers] = aeroCorridor(simInputs);
et = toc;
fprintf('COMPUTED CORRIDOR IN %2.4f s\n', et)


fitname = 'poly5';
fitlb = fit(simInputs.BC_range', Gammas(:,1), fitname);
fitub = fit(simInputs.BC_range', Gammas(:,2), fitname);
fits = {fitlb, fitub};

AllMin = min(DesignDrivers(:,1:4), [], 1);
AllMax = max(DesignDrivers(:,1:4), [], 1);

AllMin(2) = 0;
AllMax(AllMax == 0) = 1;
minmax = [AllMin',AllMax'];


simInputs.Opti.options.PlotFcn = [];
simInputs.Opti.options.Display = 'none';




%% SOO & MOO EXECUTION TIME with fits
x0 = [40.6000	-0.1765];

nIter = 10;

timeArrray = zeros([nIter, 1]);
w = [1, 1, 1, 1, 0, 0, 0;
    1, 1, 0, 0, 0, 0, 0;
    1, 0, 1, 0, 0, 0, 0;
    5, 0, 1, 0, 0, 0, 0;
    10, 5, 1, 0, 0, 0, 0];
for mode = 1:5
    simInputs.Opti.weights = w(mode,:);
    disp(simInputs.Opti.weights)


    funcCount = 0;
    for i = 1:nIter
        tic;
        [sooResults] = optimiseGD(BC_range, fits, DesignDrivers, simInputs, x0);
        timeArrray(i) = toc;
        funcCount = funcCount + sooResults.optiOutput.funcCount;
    end
    fprintf("SOO =>\tFunction counts: %3i", funcCount)
    fprintf(",   Time per iter: %10.8f", sum(timeArrray)/nIter)
    fprintf(",   Time per func. call: %10.8f\n", sum(timeArrray)/funcCount)

    fprintf("SOO optima (BC, FPA) = (%6.2f, %6.4f),\t",sooResults.BCOpt, sooResults.GammaOpt*180/pi)
    disp(sooResults.designDrivers(1:4))
    funcCount = 0;
    for i = 1:nIter
        tic;
        [mooResults] = optimiseMoo(BC_range', fits, DesignDrivers, simInputs, x0);
        timeArrray(i) = toc;
        funcCount = funcCount + mooResults.optiOutput.funcCount;
    end
    fprintf("MOO =>\tFunction counts: %3i", funcCount)
    fprintf(",   Time per iter: %10.8f", sum(timeArrray)/nIter)
    fprintf(",   Time per func. call: %10.8f\n", sum(timeArrray)/funcCount)

    fprintf("MOO optima (BC, FPA) = (%6.2f, %6.4f),\t",mooResults.BCOpt, mooResults.GammaOpt*180/pi)
    disp(mooResults.designDrivers(1:4))

    fprintf("\n")
end
disp("DONE COMPUTING SOO AND MOO WITH FITS OPTIMIZATION")
% %% SOO TOLERANCE ANALYSIS
%
% tolerances = logspace(-4,-16, 8);
% nTol = length(tolerances);
%
% optimal = zeros([nTol,2]);
% costs = zeros([nTol,1]);
% timeArrray = zeros([nTol, 1]);
%
% for i = 1:nTol
%     simInputs.Opti.options.FunctionTolerance = tolerances(i);
%     simInputs.Opti.options.OptimalityTolerance = tolerances(i);
%     tic;
%     [optivalsoo, sooResults, ~] = optimiseGD(BC_range, Gammas, DesignDrivers, simInputs, x0);
%     timeArrray(i) = toc;
%     costs(i) = sooResults.Cost;
%     optimal(i,:) = [sooResults.BCOpt, sooResults.GammaOpt*180/pi];
%
%     fprintf("SOO optima (BC, FPA) = (%6.2f, %6.4f),   Cost = %6.4f",sooResults.BCOpt, sooResults.GammaOpt*180/pi, sooResults.Cost)
%     fprintf(",   Iterations: %3i", sooResults.optiOutput.iterations)
%     fprintf(",   Function counts: %3i", sooResults.optiOutput.funcCount)
%     fprintf(",   Elpased time: %10.8f", timeArrray(i))
%     fprintf("\n")
% end

%% SOO & MOO EXECUTION TIME with NO fits
nIter = 1;

timeArrray = zeros([nIter, 1]);

for mode = 1:1
    if mode  == 1
        simInputs.Opti.weights =[ 1, 1, 1, 1, 0, 0, 0];
    else
        simInputs.Opti.weights = zeros([1,7]);
        simInputs.Opti.weights(mode) = 1;
    end
    disp(simInputs.Opti.weights)
    funcCount = 0;
    for i = 1:nIter
        tic;
        [optivalsoo, sooResults, ~] = optimiseGD_noFits(BC_range, DesignDrivers, simInputs, x0);
        timeArrray(i) = toc;
        funcCount = funcCount + sooResults.optiOutput.funcCount;
    end
    fprintf("SOO =>\tFunction counts: %3i", funcCount)
    fprintf(",   Time per iter: %10.8f", sum(timeArrray)/nIter)
    fprintf(",   Time per func. call: %10.8f\n", sum(timeArrray)/funcCount)

    funcCount = 0;
    for i = 1:nIter
        tic;
        [optivalmoo, mooResults] = optimiseMoo_noFits(BC_range', simInputs, x0);
        timeArrray(i) = toc;
        funcCount = funcCount + mooResults.optiOutput.funcCount;
    end
    fprintf("MOO =>\tFunction counts: %3i", funcCount)
    fprintf(",   Time per iter: %10.8f", sum(timeArrray)/nIter)
    fprintf(",   Time per func. call: %10.8f\n", sum(timeArrray)/funcCount)

    fprintf("\n")
end

numTol = 20;
tol = logspace(-3, -8, numTol);
nBC = 7;
BC_range = linspace(3, 60, nBC);

fitCell = cell(numTol,2);
time_array = zeros([numTol,1]);
for it = 1:numTol
    simInputs.tol = tol(it);
    Gammas = zeros([nBC,2]);
    tic;
    for i = 1:nBC
        % Compute lower and upper robust FPA angles (binary search)
        [Gammas(i,1), Gammas(i,2)] = getRobustCorridor(simInputs, BC_range(i));

        Dgamma= (Gammas(i,2) - Gammas(i,1))/2;

        % Running LOWER boundary trajectory with HIGH density mode
        % to get driver values
        [simParams, ~] = Trajectory3DSim([BC_range(i),Gammas(i,1)], 3, simInputs);
        [~, temp] = cost_function(simInputs, simParams);
        % Add simualtion results:
        DesignDrivers(i,1:7) = [temp(1:end-1), Dgamma];

        % Running UPPER boundary trajectory with LOW density mode
        % to get driver values
        [simParams, ~] = Trajectory3DSim([BC_range(i),Gammas(i,2)], 2, simInputs);

        [~, temp] = cost_function(simInputs, simParams);
        % Add simualtion results:
        DesignDrivers(i,8:14) = [temp(1:end-1), Dgamma];
    end
    time_array(it) = toc;
    fprintf('Done corridor %i of %i \n', it, numTol)

    fitname = 'poly6';
    fitCell{it,1} = fit(BC_range', Gammas(:,1), fitname);
    fitCell{it,2} = fit(BC_range', Gammas(:,2), fitname);
end
figure;
loglog(tol, time_array)
% COMPUTE EARTH CORRIDOR

BCArray = linspace(3, 60, 570);
nBC = length(BCArray);
Gammas = zeros([nBC,2]);

parfor i = 1:nBC
    % Compute lower and upper robust FPA angles (binary search)
    tempG = zeros([2,1]);
    [tempG(1), tempG(2)] = getRobustCorridor(simInputs, BCArray(i));
    Dgamma= (tempG(2) - tempG(1))/2;

    % Running LOWER boundary trajectory with HIGH density mode
    % to get driver values
    [simParams, ~] = Trajectory3DSim([BCArray(i),tempG(1)], 3, simInputs);
    [~, temp] = cost_function(simInputs, simParams);

    % Running UPPER boundary trajectory with LOW density mode
    % to get driver values
    [simParams, ~] = Trajectory3DSim([BCArray(i),tempG(2)], 2, simInputs);
    [~, temp] = cost_function(simInputs, simParams);
    Gammas(i,:) = tempG; 
    fprintf("Done %i / %i\n", i, nBC)
end

%%
close all
errorLow = zeros([length(BCArray), numTol]);
errorUp = zeros([length(BCArray), numTol]);
for i = 1:numTol
    errorLow(:,i) = (Gammas(:,1) - fitCell{i,1}(BCArray));
    errorUp(:,i) = (Gammas(:,2) - fitCell{i,2}(BCArray));
end

step = 5;
colors = colororder(parula(ceil(numTol/step)));
x = linspace(BC_range(1), BC_range(end), 100);
fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 16);
set(fig,'defaulttextinterpreter','tex');

hold on
for i = 1:step:numTol
    plot(x, 180/pi*fitCell{i,1}(x), DisplayName=sprintf("Tolerance = %1.1e", tol(i)), color=colors(ceil(i/step),:))
    plot(x, 180/pi*fitCell{i,2}(x), HandleVisibility="off", color=colors(ceil(i/step),:))
end
plot(BCArray(1:10:end), Gammas(1:10:end,1)*180/pi, '--r', DisplayName='Analytical')
plot(BCArray(1:10:end), Gammas(1:10:end,2)*180/pi, '--r', HandleVisibility='off')
legend();
grid on; grid minor
axis square;
title("Polyfits")
xlabel("\beta")
ylabel("\gamma (\circ)")


fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 16);
set(fig,'defaulttextinterpreter','tex');
title("Lower fit error")
xlabel("\beta")
ylabel("Error (%)")
axis square
grid on; grid minor
hold on;
for i = 1:step:numTol
    plot(BCArray, errorLow(:,i)./Gammas(:,1)*100, DisplayName=sprintf("Tol %1.3e", tol(i)), color=colors(ceil(i/step),:))
end
legend()


fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 16);
set(fig,'defaulttextinterpreter','tex');title("Upper fit error")
xlabel("\beta")
ylabel("Error (%)")
grid on; grid minor
axis square;
hold on;
for i = 1:step:numTol
    plot(BCArray, errorUp(:,i)./Gammas(:,2)*100, DisplayName=sprintf("Tol %1.3e", tol(i)), color=colors(ceil(i/step),:))
end

fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 12);
set(fig,'defaulttextinterpreter','tex');
yyaxis right
loglog(tol, vecnorm(errorUp,2,1)*180/pi, '-x', DisplayName='Upper Corridor')
hold on
loglog(tol, vecnorm(errorLow,2,1)*180/pi, '--o', DisplayName='Lower Corridor')
ylabel("L2 norm (\circ)")
ylim([0.5e-1,inf])
yyaxis left
loglog(tol, time_array, '.-', HandleVisibility='off')
ylabel("Execution time (s)")

grid on; grid minor;
axis square
xlabel("Binary search tolerance")
xlim([-inf,1e-3])
legend(EdgeColor='none', Color='none')









%% ACCURACY OF FIT WRT TO TOL OF BINARY SEARCH

numTol = 20;
tol = logspace(-3, -8, numTol);
nBC = 7;
BC_range = linspace(3, 60, nBC);

fitCell = cell(numTol,2);
time_array = zeros([numTol,1]);
for it = 1:numTol
    simInputs.tol = tol(it);
    Gammas = zeros([nBC,2]);
    tic;
    for i = 1:nBC
        % Compute lower and upper robust FPA angles (binary search)
        [Gammas(i,1), Gammas(i,2)] = getRobustCorridor(simInputs, BC_range(i));

        Dgamma= (Gammas(i,2) - Gammas(i,1))/2;

        % Running LOWER boundary trajectory with HIGH density mode
        % to get driver values
        [simParams, ~] = Trajectory3DSim([BC_range(i),Gammas(i,1)], 3, simInputs);
        [~, temp] = cost_function(simInputs, simParams);
        % Add simualtion results:
        DesignDrivers(i,1:7) = [temp(1:end-1), Dgamma];

        % Running UPPER boundary trajectory with LOW density mode
        % to get driver values
        [simParams, ~] = Trajectory3DSim([BC_range(i),Gammas(i,2)], 2, simInputs);

        [~, temp] = cost_function(simInputs, simParams);
        % Add simualtion results:
        DesignDrivers(i,8:14) = [temp(1:end-1), Dgamma];
    end
    time_array(it) = toc;
    fprintf('Done corridor %i of %i \n', it, numTol)

    fitname = 'poly6';
    fitCell{it,1} = fit(BC_range', Gammas(:,1), fitname);
    fitCell{it,2} = fit(BC_range', Gammas(:,2), fitname);
end
figure;
loglog(tol, time_array)
%%



load("C:\GIT\D-AsTRO-research\Output files\savedExperiment_testEarthVinf_1p0.mat")
close all
errorLow = zeros([length(BCArray), numTol]);
errorUp = zeros([length(BCArray), numTol]);
for i = 1:numTol
    errorLow(:,i) = (Gammas(:,5) - fitCell{i,1}(BCArray));
    errorUp(:,i) = (Gammas(:,4) - fitCell{i,2}(BCArray));
end

step = 5;
colors = colororder(parula(ceil(numTol/step)));
x = linspace(BC_range(1), BC_range(end), 100);
fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 16);
set(fig,'defaulttextinterpreter','tex');

hold on
for i = 1:step:numTol
    plot(x, 180/pi*fitCell{i,1}(x), DisplayName=sprintf("Tolerance = %1.1e", tol(i)), color=colors(ceil(i/step),:))
    plot(x, 180/pi*fitCell{i,2}(x), HandleVisibility="off", color=colors(ceil(i/step),:))
end
plot(BCArray(1:10:end), Gammas(1:10:end,5)*180/pi, '--r', DisplayName='Analytical')
plot(BCArray(1:10:end), Gammas(1:10:end,4)*180/pi, '--r', HandleVisibility='off')
legend();
grid on; grid minor
axis square;
title("Polyfits")
xlabel("\beta")
ylabel("\gamma (\circ)")


fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 16);
set(fig,'defaulttextinterpreter','tex');
title("Lower fit error")
xlabel("\beta")
ylabel("Error (%)")
axis square
grid on; grid minor
hold on;
for i = 1:step:numTol
    plot(BCArray, errorLow(:,i)./Gammas(:,5)*100, DisplayName=sprintf("Tol %1.3e", tol(i)), color=colors(ceil(i/step),:))
end
legend()
ylim([-0.4, 1])


fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 16);
set(fig,'defaulttextinterpreter','tex');title("Upper fit error")
xlabel("\beta")
ylabel("Error (%)")
grid on; grid minor
axis square;
hold on;
for i = 1:step:numTol
    plot(BCArray, errorUp(:,i)./Gammas(:,4)*100, DisplayName=sprintf("Tol %1.3e", tol(i)), color=colors(ceil(i/step),:))
end
ylim([-0.4, 1])

fig = figure();
set(fig, 'defaultlinelinewidth',2);
set(fig, 'defaultaxesfontsize', 12);
set(fig,'defaulttextinterpreter','tex');
yyaxis right
loglog(tol, vecnorm(errorUp,2,1)*180/pi, '-x', DisplayName='Upper Corridor')
hold on
loglog(tol, vecnorm(errorLow,2,1)*180/pi, '--o', DisplayName='Lower Corridor')
ylabel("L2 norm (\circ)")

yyaxis left
loglog(tol, time_array, '.-', HandleVisibility='off')
ylabel("Execution time (s)")

grid on; grid minor;
axis square
xlabel("Binary search tolerance")
xlim([-inf,1e-3])
legend()
