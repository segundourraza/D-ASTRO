function [y, designDrivers, orbitalElements] = MultiObejctiveFunc(x, simInputs, minmax)

density_mode = 3;


[simParams, ~] = Trajectory3DSim(x, density_mode, simInputs);
orbitalElements = simParams([1,2,3,9, 10, 11]);

if nargin > 2
    [~, designDrivers] = cost_function(simInputs, simParams, minmax);
else
    [~, designDrivers] = cost_function(simInputs, simParams);
end
y = designDrivers(1:4);
% fprintf("%f,   %f\n", x(1), x(2))
% fprintf("%f,   %f,   %f,   %f\n", y(1), y(2), y(3), y(4))
% y(2) = y(2)/1e3;
% y(3) = y(3)/1e6;
% y(4) = y(4)/1e4;