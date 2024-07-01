clc
clc
clear
close all

cd("..\")

plotOption = false;
saveOption = false;

AerocaptureCorridor = 2;
densityMode = 1;

simInputs = inputSim();
BC_range = linspace(3, 60, 7);

optimisation = 1;
weights = [10, 1, 1, 1]';
targetOrbit = [ 4620.5;70;0.05];
run caseSetup.m

[Gammas, OrbitalParams, Trajectories, DesignDrivers] = aeroCorridor(simInputs);


nTheta = 7;
mu = zeros([1,nTheta]);
sigma = zeros([1,nTheta]);
AllMax = zeros([1,nTheta]);
AllMin = zeros([1,nTheta]);

for i = 1:nTheta
    % Standerisation
    mu(i) = mean(DesignDrivers(:,i:nTheta:end), 'all');
    sigma(i) = std(DesignDrivers(:,i:nTheta:end), [], 'all');

    % Normalisation
    AllMin(i) = min(real(DesignDrivers(:,i:nTheta:end)), [],'all');
    AllMax(i) = max(real(DesignDrivers(:,i:nTheta:end)), [],'all');
end

simInputs.Opti.mins = AllMin;
simInputs.Opti.mins(2) = 0;
simInputs.Opti.maxs = AllMax;


cd("Performance\")

%% GERNERATING ALLCOST MESH (COMPUTATIONALLY EXPENSIEVE)
BCArray = linspace(3, 60, 570);

nfpa = 3e2;
maxg = max(Gammas(:,2));
ming = min(Gammas(:,1));

fpa = linspace(ming, maxg, nfpa);
[BCMesh, FPA] = meshgrid(BCArray, fpa);

allCost = zeros(size(BCMesh));
allDV =  zeros(size(BCMesh));
allQ = zeros(size(BCMesh));
allqmax = zeros(size(BCMesh)); 

allVelAtHmin = zeros(size(BCMesh)); 

tracker = 0;
parfor ibc = 1: length(BCArray)
% for ibc = 1: length(BCArray)
    fpa = linspace(ming, maxg, nfpa);
    temp = zeros([1,nfpa]);
    dv = zeros([1, nfpa]);
    Q = zeros([1, nfpa]);
    qmax = zeros([1, nfpa]);
    orbitalElements = zeros([1. nfpa]);
    for ifpa = 1:nfpa
        [temp(ifpa), drivers, coe] = ObejctiveFunc2D([BCArray(ibc), fpa(ifpa)],simInputs);
        dv(ifpa) = drivers(2);
        Q(ifpa) = drivers(3);
        qmax(ifpa) = drivers(4);

        % temp(ifpa) = sqrt(BCArray(ibc)/cos(-fpa(ifpa)));
        % Q(ifpa) = temp(ifpa);

        % temp(ifpa) = sqrt(BCArray(ibc)*sin(-fpa(ifpa)));
    end
    allCost(:,ibc) = temp;
    allDV(:,ibc) = dv;
    allQ(:,ibc) = Q;
    allqmax(:,ibc) = qmax;
    fprintf("DONE " + ibc+ " / "  + length(BCArray) +  ".\n")
end

save('volume10_allCost.mat','allCost', "BCArray", "fpa")
save('volume10_allDV.mat','allDV', "BCArray", "fpa")
save('volume10_allQ.mat','allQ', "BCArray", "fpa")
save('volume10_allqmax.mat','allqmax', "BCArray", "fpa")