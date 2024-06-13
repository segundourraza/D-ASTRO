function [cost] = least_square_method(x, y, k)

cost = sum((k.*(x-y)).^2,2);
end