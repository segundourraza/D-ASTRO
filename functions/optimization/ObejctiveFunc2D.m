function [y,designDrivers, orbitalElements] = ObejctiveFunc2D(x, simInputs)

density_mode = 3;

[simParams, ~] = Trajectory3DSim(x, density_mode, simInputs);
orbitalElements = simParams([1,2,3,9, 10, 11]);

[y, designDrivers] = cost_function(simInputs, simParams, 1);
% [y, designDrivers] = cost_function(simInputs, simParams);
