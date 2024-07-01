clc
clear
close all

cd("..\")

plotOption = false;
saveOption = false;

AerocaptureCorridor = 2;
densityMode = 1;

BC_range = linspace(3, 60, 7);

optimisation = 1;
weights = [0, 0, 1, 0]';
targetOrbit = [ 4620.5;70;0.05];
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

%% SOO TOLERANCE ANALYSIS
clc


x0 = [59.0000	-1.0472;
    40.6000	-0.1765;
    59.0000	-0.0010;
    3.0010	-1e-3;
    39.6000	-0.1894;
    37.8000	-0.1833];
%
% x0 = x0(3,:);
nx0 = size(x0, 1);

mintol = -10; % default step tol

tolerances = logspace(mintol,-16, 8);
tolerances = logspace(mintol,-16, 4);

nTol = length(tolerances);

allXf = cell(1,nx0);
% parfor j = 1:size(x0,1)
for j = 1:size(x0,1)
    localsimInputs = simInputs;
    fprintf("X0: (%6.2f, %6.4f)\n",x0(j,1), x0(j,2)*180/pi)
    Xf = zeros([nTol, 2]);
    for i = 1:nTol
        localsimInputs.Opti.options.FunctionTolerance = tolerances(i);
        localsimInputs.Opti.options.StepTolerance = tolerances(i);
        localsimInputs.Opti.options.OptimalityTolerance = tolerances(i);
        [sooResults] = optimiseGD(BC_range, fits, DesignDrivers, localsimInputs, x0(j,:));
        Xf(i,:) = [sooResults.BCOpt, sooResults.GammaOpt];
        fprintf("Step Tol = %1.2e", localsimInputs.Opti.options.StepTolerance)
        fprintf(",   Func Tol = %1.2e", localsimInputs.Opti.options.FunctionTolerance)
        fprintf(",   (%6.2f, %6.4f)",sooResults.BCOpt, sooResults.GammaOpt*180/pi)
        fprintf(",   Q = %2.3e", sooResults.designDrivers(3))
        fprintf(",   Iterations: %3i", sooResults.optiOutput.iterations)
        fprintf(",   Function counts: %3i", sooResults.optiOutput.funcCount)
        fprintf("\n")
    end
    allXf{j}= Xf;
end
disp("DONE")

close all
fig = figure();
x = linspace(3, 60, 100);
hold on
plot(x, fits{1}(x)*180/pi, '-m')
plot(x, fits{2}(x)*180/pi, '-m')
yl = ylim();

ls = {'-.', '--',':','-' };

for i = 1:nx0
    for j = 4
        text(x0(i,1), x0(i,2)*180/pi, string(i), FontSize=24)
        plot([x0(i,1), allXf{i}(j,1)], [x0(i,2),allXf{i}(j,2)]*180/pi, LineStyle=ls{j})
    end
end
ylim(yl);
%%

BC = 3;
fpa =linspace(fits{1}(BC), fits{2}(BC), 100);

simInputs3 = simInputs;
simInputs3.AerocaptureCorridor = 1;
simInputs3.densityMode = 3;
simInputs3.BC_range = BC*ones([1,length(fpa)]);
simInputs3.gamma_range = fpa;

[Gammas, OrbitalParams, Trajectories, DesignDrivers2] = aeroCorridor(simInputs3);

figure;
hold on
yyaxis left
plot(fpa*180/pi, DesignDrivers2(:,3)/1e6)
xline(fits{1}(BC)*180/pi + simInputs3.Dgamma*180/pi)
xline(fits{2}(BC)*180/pi - simInputs3.Dgamma*180/pi)

for i = 1:nx0
    xline(allXf{i}(end,2)*180/pi, 'g')
end


yyaxis right
plot(fpa*180/pi, (DesignDrivers2(:,3)-AllMin(3))/(AllMax(3) - AllMin(3)))








%% SOO & MOO EXECUTION TIME with fits

x0 = [40.6000	-0.1765];
x0 = [39.6000	-0.1894];

nIter = 10;
scaling = [ 1.15, 5.15;
    0, 3.56;
    3.42, 41.53;
    3.45, 46.32];
scaling = scaling';

timeArrray = zeros([nIter, 1]);
w = [1, 1, 1, 1, 0, 0, 0;
    1, 1, 0, 0, 0, 0, 0;
    1, 0, 1, 0, 0, 0, 0;
    5, 0, 1, 0, 0, 0, 0;
    0.1, 0, 100, 0,0,0, 0;
    10, 5, 1, 0, 0, 0, 0];
for mode = 1:6
    simInputs.Opti.weights = w(mode,:);
    disp(simInputs.Opti.weights(1:4))


    funcCount = 0;
    for i = 1:nIter
        tic;
        [sooResults] = optimiseGD(BC_range, fits, DesignDrivers, simInputs, x0);
        timeArrray(i) = toc;
        funcCount = funcCount + sooResults.optiOutput.funcCount;
    end
    fprintf("SOO =>\tTime per iter: %10.8f", sum(timeArrray)/nIter)
    fprintf(",   Function counts: %3i", funcCount)
    fprintf(",   Time per func. call: %10.8f ms\n", sum(timeArrray)/funcCount*1e3)

    fprintf("SOO optima (BC, FPA) = (%6.2f, %6.4f),\n",sooResults.BCOpt, sooResults.GammaOpt*180/pi)
    fprintf(['Drivers:      ' repmat(' %2.4f ',1,numel(sooResults.designDrivers(1:4))) '\n'],sooResults.designDrivers(1:4));
    d = (sooResults.designDrivers(1:4) - scaling(1,:))./(scaling(2,:) - scaling(1,:));
    fprintf(['Norm Drivers: ' repmat(' %2.4f ',1,numel(d)) '\n\n'],d);


    funcCount = 0;
    for i = 1:nIter
        tic;
        [mooResults] = optimiseMoo(BC_range, fits, DesignDrivers, simInputs, x0);
        timeArrray(i) = toc;
        funcCount = funcCount + mooResults.optiOutput.funcCount;
    end
    fprintf("MOO =>\tTime per iter: %10.8f", sum(timeArrray)/nIter)
    fprintf(",   Function counts: %3i", funcCount)
    fprintf(",   Time per func. call: %10.8f ms\n", sum(timeArrray)/funcCount*1e3)

    fprintf("MOO optima (BC, FPA) = (%6.2f, %6.4f),\n",mooResults.BCOpt, mooResults.GammaOpt*180/pi)
    fprintf(['Drivers:      ' repmat(' %2.4f ',1,numel(mooResults.designDrivers(1:4))) '\n'],mooResults.designDrivers(1:4));
    d = (mooResults.designDrivers(1:4) - scaling(1,:))./(scaling(2,:) - scaling(1,:));
    fprintf(['Norm Drivers: ' repmat(' %2.4f ',1,numel(d)) '\n'],d);

    fprintf("\n")
end
disp("DONE COMPUTING SOO AND MOO WITH FITS OPTIMIZATION")







%% SOO & MOO EXECUTION TIME with NO fits



x0 = [40.6000	-0.1765];


simInputs.Opti.options.OptimalityTolerance = 1e-16;

nIter = 1;

timeArrray = zeros([nIter, 1]);

timeArrray = zeros([nIter, 1]);
w = [1, 1, 1, 1, 0, 0, 0;
    1, 1, 0, 0, 0, 0, 0;
    1, 0, 1, 0, 0, 0, 0;
    5, 0, 1, 0, 0, 0, 0;
    0.1, 0, 100, 0,0,0, 0;
    10, 5, 1, 0, 0, 0, 0];

for mode = 1:6
    simInputs.Opti.weights = w(mode,:);
    disp(simInputs.Opti.weights(1:4))

    % funcCount = 0;
    % for i = 1:nIter
    %     tic;
    %     [sooResults] = optimiseGD_noFits(BC_range, DesignDrivers, simInputs, x0);
    %     timeArrray(i) = toc;
    %     funcCount = funcCount + sooResults.optiOutput.funcCount;
    % end
    % fprintf("SOO =>\tFunction counts: %3i", funcCount)
    % fprintf(",   Time per iter: %10.8f", sum(timeArrray)/nIter)
    % fprintf(",   Time per func. call: %10.8f\n", sum(timeArrray)/funcCount)
    % fprintf("(BC, FPA) = (%6.2f, %6.4f),\n",sooResults.BCOpt, sooResults.GammaOpt*180/pi)
    % fprintf(['Drivers:      ' repmat(' %2.4f ',1,numel(sooResults.designDrivers(1:4))) '\n'],sooResults.designDrivers(1:4));

    funcCount = 0;
    for i = 1:nIter
        tic;
        [mooResults] = optimiseMoo_noFits(BC_range', simInputs, x0);
        timeArrray(i) = toc;
        funcCount = funcCount + mooResults.optiOutput.funcCount;
    end
    fprintf("MOO =>\tFunction counts: %3i", funcCount)
    fprintf(",   Time per iter: %10.8f", sum(timeArrray)/nIter)
    fprintf(",   Time per func. call: %10.8f\n", sum(timeArrray)/funcCount)
    fprintf("(BC, FPA) = (%6.2f, %6.4f),\n",mooResults.BCOpt, mooResults.GammaOpt*180/pi)
    fprintf(['Drivers:      ' repmat(' %2.4f ',1,numel(mooResults.designDrivers(1:4))) '\n'],mooResults.designDrivers(1:4));
    
    fprintf("\n")
end

%% TRAJECTORY ANALYSER

BC_range = 20;
gamma_range = linspace(-11.5, -8.5, 50)*pi/180;
k = length(gamma_range);
numTol = 10;

time_array = zeros([k,1]);
count_array = zeros([k,1]);

fig = figure();
subplot(1,2,1)
hold on;
for i = 1:k
    tic;
    for j =1 :numTol
        [simParams, trajResults] = Trajectory3DSim([BC_range,gamma_range(i)], 3, simInputs);
    end
    time_array(i) = toc/numTol;
    count_array(i) = length(trajResults(:,1));
    if mod(i,2) == 0
        plot(trajResults(:,1), (trajResults(:,2) - simInputs.R)/1e3)

        % plot(trajResults(:,1), trajResults(:,end-1)) % qmax
        fprintf("Done %3i / %3i\n",i, k)
    end
    % fprintf("%6.4f\n",gamma_range(i)*180/pi)
end
ylabel("h (km)")
xlabel("Time (s)")

subplot(1,2,2)
yyaxis left
plot(time_array*1e3)
ylabel("Execution time (ms)")

yyaxis right
plot(count_array)
ylabel("Trajectory points")

fig = figure();
plot(time_array*1e3./count_array)
ylabel("Execution time per trajectory point")

%% BINARY SEARCH ANALYSIS





numTol = 20;
tol1 = logspace(-3, -8, numTol);
BC_range = 20;

% tol1 = [2.6e-5, 1e-8];
% numTol = length(tol1);


Gamma_test = zeros([numTol,2]);
time_array = zeros([numTol,1]);
count_array = zeros([numTol,1]);

numIter = 1;
valCell = cell(1,2);
for i = 1:numTol
    simInputs.tol = tol1(i);
    tic;
    for j = 1:numIter
        BC = BC_range;
        Rplanet = simInputs.R;
        h_AI = simInputs.h_AI;
        mu = simInputs.mu;

        % START BINRAY SEARCH FOR LOWER ROBUST AEROCAPTURE CORRIDOR BOUNDARY
        % Corresponds to HIGH density mode
        params = [BC, 3, simInputs.R, simInputs.SC.m, simInputs.mu, simInputs.SC.L_D, simInputs.Omega];

        gamma_min = simInputs.gamma_min;
        gamma_max = simInputs.gamma_max;
        gamma_lower = 3*simInputs.tol;
        gamma_prev = 0;
        iteration = 1;

        val = nan([100,1]);
        while abs(gamma_lower-gamma_prev) > simInputs.tol && (iteration < simInputs.MaxIterations)
            gamma_prev = gamma_lower;
            gamma_lower = 0.5*(gamma_max + gamma_min);
            val(iteration) = gamma_lower;
            x0 = [simInputs.delta0, simInputs.lambda0, simInputs.R0, simInputs.V0, gamma_lower, simInputs.chi0];
            [~, y] = ode45( @(t, x) reentryModel(t, x, params, simInputs.atmoModel), simInputs.tspan, x0, simInputs.odeOptions);
            rEnd = y(end,3); % Radius of SC at end of trajectory

            % Identify lower limit of corridor, separates surface collision with no surface collision
            if (rEnd - Rplanet) < h_AI
                gamma_min = gamma_lower;
            else
                gamma_max = gamma_lower;
            end
            iteration = iteration + 1;
        end


        count_array(i) = iteration;
        valCell{1} = val;

        % START BINRAY SEARCH FOR UPPER ROBUST AEROCAPTURE CORRIDOR BOUNDARY
        % Corresponds to LOW density mode
        params(2) = 2;
        gamma_min = gamma_lower;
        gamma_max = simInputs.gamma_max;
        gamma_upper = 3*simInputs.tol;
        gamma_prev = 0;
        iteration = 1;
        val = nan([100,1]);
        while abs(gamma_upper-gamma_prev) > simInputs.tol && (iteration < simInputs.MaxIterations)
            gamma_prev = gamma_upper;
            gamma_upper = 0.5*(gamma_max + gamma_min);
            val(iteration) = gamma_upper;
            x0(5) =  gamma_upper;
            [~, y] = ode45( @(t, x) reentryModel(t, x, params, simInputs.atmoModel), simInputs.tspan, x0, simInputs.odeOptions);
            rEnd = y(end,3); % Altitude of SC at end of trajectory
            vEnd = y(end,4); % Velocity of SC at end of trajectory
            specific_energy = 0.5*vEnd^2 - mu/rEnd;
            % Identify lower limit of corridor, separated by negative orbital
            % specific
            if specific_energy < 0
                gamma_min = gamma_upper;
            else
                gamma_max = gamma_upper;
            end
            iteration = iteration + 1;
        end

        count_array(i) = count_array(i) +  iteration;
        valCell{2} = val;
        Gamma_test(i,:) = [gamma_lower,     gamma_upper];
    end
    fprintf("Done %i / %i\n", i, numTol)
    time_array(i) = toc/numIter;
end
a = valCell{1}(~isnan(valCell{1}));
val = [a;valCell{2}(~isnan(valCell{2}))];

fig = figure();
set(fig,'defaultlinelinewidth',2)
hold on
yline(gamma_lower*180/pi, '','Lower',LineWidth=2, Color='r',LineStyle='--')
yline(gamma_upper*180/pi, '','Upper',LineWidth=2, Color='b',LineStyle='--')
plot(val*180/pi, '-xk')
grid on; grid minor;
ylim([-14, -4])

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



load("C:\GIT\D-AsTRO-research\Output files\savedExperiment_step01V3.mat")
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
