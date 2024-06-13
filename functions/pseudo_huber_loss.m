function [cost] = pseudo_huber_loss(x, y, weight, delta)

if nargin < 4
    delta = 0.1;
end

val = (y - x);
cost = delta.^2.*(sqrt(1+sum(weight.*(val./delta).^2,2)) -1);
end