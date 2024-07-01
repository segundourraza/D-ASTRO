function [c, ceq] = binarySearchConstraint(x, Dg, simInputs)
[flo, fup] = getRobustCorridor(simInputs, x(1));
c(1) = x(2) - fup + Dg;
c(2) = flo - x(2) + Dg;
ceq = [];
end