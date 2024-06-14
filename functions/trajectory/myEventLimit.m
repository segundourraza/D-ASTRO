function [value, isterminal, direction] = myEventLimit(t, y, params)

Rplanet = params(1);
h_AI = params(2);

value      = ((y(3) >= (Rplanet+h_AI)+1) | (y(3) <= (Rplanet)));
isterminal = 1;   % Stop the integration
direction  = 0;

end