function [f] = tauberSuttonFunction(V)
% Gives the value of the tabulated function f(V) used in the Tauber-Sutton
% correlation for radiative heating

% polynomial coefficients
a0 = 1.849710801386854e+03;
a1 = -1.007772515250308;
a2 = 2.027860638323226e-04;
a3 = -1.793168028714326e-08;
a4 = 5.943163921781362e-13;


f = a0 + V.* (a1 + V.* (a2 + V.* (a3 + a4 .* V)));
f(V<6000) = 0.0;
f(V>9000) = 32.8;

end