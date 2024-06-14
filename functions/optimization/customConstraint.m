function [c, ceq] = customConstraint(x,  fup, flo, Dg)
c(1) = x(2) - fup(x(1)) + Dg;
c(2) = flo(x(1)) - x(2) + Dg;
ceq = [];
end