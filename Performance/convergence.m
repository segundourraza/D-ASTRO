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
% weights = [1, 1, 1, 1]';
weights = [0, 0, 1, 0]';
targetOrbit = [ 4620.5;70;0.05];
run caseSetup.m

cd("Performance\")
addpath("Multi-objective optimization\")

tic

% OG VERSION
[Gammas, OrbitalParams, Trajectories, DesignDrivers] = aeroCorridor(simInputs);
et = toc;
fprintf('DONE COMPUTING LIMITS IN %2.6f s\n',et)

fitname = 'poly5';
fitlb = fit(simInputs.BC_range', Gammas(:,1), fitname);
fitub = fit(simInputs.BC_range', Gammas(:,2), fitname);
fits = {fitlb, fitub};

nfpa = 3e2;
maxg = max(Gammas, [],'all');
ming = min(Gammas, [],'all');

simInputs.Opti.options.Display = 'none';
simInputs.Opti.options.PlotFcn = [];

%% CONVERGENCE PROEPRTIES
close all


x0 = [59.0000	-1.0472;
      40.6000	-0.1765;
      59.0000	-0.0010;
      3.0010	-1e-3;
      39.6000	-0.1894;
      37.8000	-0.1833];

Xf = cell(size(x0,1), 2);
for i = 1:size(x0,1)
    [sooResults] = optimiseGD(BC_range, fits, DesignDrivers, simInputs, x0(i,:));
    Xf{i,1} = [sooResults.BCOpt, sooResults.GammaOpt]; 
    fprintf("SOO (BC, FPA) = (%6.2f, %6.4f),\n",sooResults.BCOpt, sooResults.GammaOpt*180/pi)
    Xf{i,2} = [0, 0];
    % [mooResults] = optimiseMoo(BC_range, fits, DesignDrivers, simInputs, x0(i,:));
    % Xf{i,2} = [mooResults.BCOpt, mooResults.GammaOpt];
    % fprintf("MOO (BC, FPA) = (%6.2f, %6.4f),\n",mooResults.BCOpt, mooResults.GammaOpt*180/pi)

    fprintf("DONE %i of %i\n", i, size(x0,1))
end

%%
lc = 'w';
k = 14;
% load("C:\GIT\D-AsTRO-research\Output files\savedExperiment_step01V3.mat")
% load("C:\GIT\D-AsTRO-research\new variables\overall_allCost.mat")
% allCost = allCost - 0.1;

% load("volume10_allCost.mat")
load("results_allCost.mat")
% load("simualtionParams.mat")

BCArray = linspace(3, 60, 570);
fpa = linspace(ming, maxg, nfpa);
[BCMesh, FPA] = meshgrid(BCArray, fpa);

fig =figure();
set(fig, 'defaultaxesfontsize', k)
set(fig, 'defaultlinelinewidth', 2.2)
colormap(hot)

% [C, h] = contourf(BCMesh,FPA.*180/pi,allCost, [0,1.0,1.2,1.3, 1.4, 2, 4, 6, 8, 10, 20, 40 ],'EdgeColor',lc);
[C, h] = contourf(BCMesh,FPA.*180/pi,allCost,[0, 0.4, 0.5,0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8], 'EdgeColor',lc);
set(gca, 'ColorScale', 'log')

hold on
plot(BCArray, fits{1}(BCArray)*180/pi, '-m')
plot(BCArray, fits{2}(BCArray)*180/pi, '-m')

for i = 1:size(x0,1)
    plot([x0(i,1), Xf{i,1}(1)], [x0(i,2), Xf{i,1}(2)]*180/pi, '-g')
    plot(Xf{i,1}(1), Xf{i,1}(2)*180/pi, 'xg')
    plot(x0(i,1), x0(i,2)*180/pi, 'ko', MarkerFaceColor='k')
end


for i = 1:size(x0,1)
    plot([x0(i,1), Xf{i,2}(1)], [x0(i,2), Xf{i,2}(2)]*180/pi, '-b')
    plot(Xf{i,2}(1), Xf{i,2}(2)*180/pi, 'xb')
end

ylim([ming, maxg]*180/pi)