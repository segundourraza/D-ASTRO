function [f] = tauberSuttonFunction(V)
% Gives the value of the tabulated function f(V) used in the Tauber-Sutton
% correlation for radiative heating

% polynomial coefficients
p1 = 1.808638254121861e-19;
p2 = -7.749778685483423e-15;
p3 = 1.380022203185181e-10;
p4 = -1.308026283879384e-06;
p5 = 0.006965644994364;
p6 = -19.770449216103156;
p7 = 2.336639165179594e+04;

if V < 6000
    f = 0.0;
elseif V >9000
    f = 32.8;
else
    f = p7 + V.* (p6 + V.* (p5 + V.* (p4 +  V.*(p3 + V.*(p2 + V.*p1)))));
end
end