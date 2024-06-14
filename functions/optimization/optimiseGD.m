function [Optimalvalues, Results,simInputs] = optimiseGD(BC_Array, Gammas, DesignDrivers, simInputs, x0)
if nargin <5
    ind = floor(length(BC_Array)/2);
    x0 = [BC_Array(ind), (Gammas(ind,2) + Gammas(ind,1))/2];
    x0 = [BC_Array(ind), Gammas(ind,2)];
end
simInputs.BC_range = BC_Array;
Gammalb = Gammas(:,1);
Gammaub = Gammas(:,2);

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

fitname = 'poly5';
fitlb = fit(BC_Array', Gammalb, fitname);
fitub = fit(BC_Array', Gammaub, fitname);
simInputs.fitub = fitub;
simInputs.fitlb = fitlb;
%% OPTIMISATION

options = simInputs.Opti.options;

fun = @(x) ObejctiveFunc2D(x, simInputs);
nonlcon = @(x) customConstraint(x, fitub, fitlb, simInputs.Dgamma);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [BC_Array(1), -pi/2];
ub = [BC_Array(end), -0];
xOpt = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
[cost,designDrivers, COE] = ObejctiveFunc2D(xOpt, simInputs);

designDrivers(2) = designDrivers(2)/1e3; % m/s to Km/s
designDrivers(3) = designDrivers(3)/1e6; % J/m2 to MJ/m2
designDrivers(4) = designDrivers(4)/1e4; % W/m2 to W/cm2;
designDrivers(6) = designDrivers(6)/1e3; % m to km

C = [cost, xOpt(1), xOpt(2), designDrivers];
COE_target = [simInputs.Opti.a_orb_target,  simInputs.Opti.e_orb_target, simInputs.Opti.i_orb_target, nan, nan, nan];

%% TABLES
if simInputs.verbose
    var_names = {'Cost','BC','Flight path angle','Volume umbrella [m3]', 'Delta V[Km/s]','Heat Load[MJ/cm2]', 'Max Heat transfer[W/cm2]', 'Max Wall Temp.[K]', 'Delta Gamma [deg]' };
    T1 = table(C(:,1),C(:,2), C(:,3)*180/pi, C(:,4), C(:,5), (C(:,6))*180/pi, C(:,7), C(:,8), C(:,10)*180/pi, 'VariableNames', var_names);
    disp(T1)
    
    var_names = {'Semi-Major axis[Km]', 'eccenreicity', 'Inlcination', 'Argument of Right Ascending node', 'Argument of perigee', 'True Anomaly'};
    row_names = {'Actual', 'target'};
    D2 = [COE;COE_target];
    D2(:,1) = D2(:,1)/1e3;
    T2 = table(D2(:,1), D2(:,2), D2(:,3)*180/pi, D2(:,4)*180/pi, (D2(:,5))*180/pi, D2(:,6)*180/pi, ...
        'VariableNames', var_names, ...
        'RowNames',row_names);
    T2 = varfun(@(x) num2str(x, ['%' sprintf('.%df', 2)]), T2);
    T2.Properties.VariableNames = var_names;
    T2.Properties.RowNames = row_names;
    disp(T2)
end
%% OUTTPUTS
Optimalvalues = [C, COE];
Results.Cost = C(1);
Results.BCOpt = C(2);
Results.GammaOpt = C(3);
Results.designDrivers = designDrivers;
Results.COE = COE;
Results.x0 = x0;
Results.target = simInputs.Opti.target;
Results.Weights = simInputs.Opti.weights;
Results.mins = simInputs.Opti.mins;
Results.maxs = simInputs.Opti.maxs;
Results.fits = {fitlb, fitub};

end