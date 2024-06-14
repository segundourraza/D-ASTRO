function [a] = axial_acceleration(x, BC, simInputs)    

m = simInputs.SC.m;
mu = simInputs.mu;
Rplanet = simInputs.R;
L_D = simInputs.SC.L_D;


density_mode = 3;
Omega = 2*pi / (24*3600 + 38*60 + (22+35)/2);

R = x(:,1);
V = x(:,2);
gamma = x(:,3);
lambda = x(:,4); %latitude
chi = x(:,5);
A = pi/2 - chi ;

% Gravitational accelertaion 
g = mu./R.^2;

% Atmopsheric parameters

if R < Rplanet
    R = Rplanet;
elseif R > Rplanet + 200e3
    R = Rplanet +200e3;
end

h = R - Rplanet; % Altitude

[~, ~, rho] = marsatm(h, density_mode, 0, simInputs.models); % Atmospheric properties
q_inf = 0.5.*rho.*V.^2; % Dynamic head

%% DRAG MODEL
D = q_inf .* m./BC ;

%% LIFT MODEL
L = D.*L_D;

%% FORCES
Xfo = L;
Zfo = 0;

%% AXIAL ACCELERATIONS

phi_dot = (V./R-g./V).*cos(gamma) + Xfo./(V.*m) + 2*Omega.*cos(A).*cos(lambda)+ ...
    (Omega^2.*R.*cos(lambda).*(cos(gamma).*cos(lambda)+sin(gamma).*sin(A).*sin(lambda)))./V;
heading_dot = V.*sin(A).*tan(lambda).*cos(gamma)./R-g.*sin(A)./V-Zfo./(V.*cos(gamma).*m)-...
    2.*Omega.*(tan(gamma).*cos(A).*cos(lambda)-sin(lambda))+(Omega.^2.*R.*sin(A).*sin(lambda).*cos(lambda))./(V.*cos(gamma));

a = V.*sqrt(phi_dot.^2+heading_dot.^2);
