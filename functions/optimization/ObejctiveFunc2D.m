function [y,designDrivers, orbitalElements] = ObejctiveFunc2D(x, simInputs, minmax)

density_mode = 3;

[simParams, ~] = Trajectory3DSim(x, density_mode, simInputs);
orbitalElements = simParams([1, 2, 3, 9, 10, 11]);

if nargin > 2
    [y, designDrivers] = cost_function(simInputs, simParams, minmax);
else
    [y, designDrivers] = cost_function(simInputs, simParams);
end
