function [f] = tauberSuttonFunction_earth(V)
% Gives the value of the tabulated function f(V) used in the Tauber-Sutton
% correlation for radiative heating

% polynomial coefficients
a4 = -7.455362759059251e-13;
a3 = 3.563683167923780e-08;
a2 = -5.857080241069495e-04;
a1 = 4.073786910844869;
a0 = -1.031096974052802e+04;

f = a0 + V.* (a1 + V.* (a2 + V.* (a3 + a4 .* V)));
f(V<9000) = 0.0;
f(V>16e3) = 2040;

end