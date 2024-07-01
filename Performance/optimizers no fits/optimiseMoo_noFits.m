function [Results] = optimiseMoo_noFits(BC_Array, simInputs, x0)
simInputs.BC_range = BC_Array;

%% OPTIMISATION


options = optimoptions('fgoalattain', ...
                       'display', 'none',...
                       'MaxIterations',3500, ...
                       'FunctionTolerance',1e-16, 'StepTolerance',1e-16);
options.OptimalityTolerance = 1e-16;
fun = @(x) MultiObejctiveFunc(x, simInputs);
nonlcon = @(x) binarySearchConstraint(x, simInputs.Dgamma, simInputs);


weight = simInputs.Opti.weights(1:4);
weight(weight ~= 0) = 1./weight(weight ~= 0);
goal = [0, 0, 0, 0];

A = [];
b = [];
Aeq = [];
beq = [];
lb = [BC_Array(1), -pi/2];
ub = [BC_Array(end), -0];
[xOpt,~,~,~,output,~] = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlcon,options);
[cost,designDrivers, COE] = MultiObejctiveFunc(xOpt, simInputs);

designDrivers(2) = designDrivers(2)/1e3; % m/s to Km/s
designDrivers(3) = designDrivers(3)/1e6; %  J/m2 to MJ/cm2
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
    row_names = {'Actual', 'Target'};
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
Results.Cost = C(1);
Results.BCOpt = xOpt(1);
Results.GammaOpt = xOpt(2);
Results.designDrivers = designDrivers;
Results.COE = COE;
Results.x0 = x0;
Results.Target = simInputs.Opti.target;
Results.Weights = simInputs.Opti.weights;
Results.optiOutput = output;
end