function [y] = reentryModel(t, x, params, models)

%% GLOBAL PARAMETERS
BC = params(1);
density_mode = params(2);
Rm = params(3);
m = params(4);
mu = params(5);
L_D = params(6);
Omega = params(7);

delta = x(1); % Longitude
lambda = x(2); % Latitude
R = x(3); % Position 
V = x(4); % Velocity
gamma = x(5); % FPA
chi = x(6); % HEading angle (pi/2 - A)

w3 = 0;

% Gravitational accelertaion 
g = mu/R^2;

% Atmopsheric parameters

if R < Rm
    R = Rm;
elseif R > Rm + 200e3
    R = Rm +200e3;
end

h = R - Rm; % Altitude
[~, ~, rho] = atmoModel(h, density_mode, 0, models); % Atmospheric properties
% [~, ~, rho] = earthAtm(h, density_mode);
q_inf = 0.5*rho*V^2; % Dynamic head

%% DRAG MODEL
D = q_inf .* m./BC;

%% LIFT MODEL
L = D.*L_D;

%% FORCES
Xfo = L;
Yfo = -D;
Zfo = 0;

%% GOVERNING EQUATIONS

y(1) = V*cos(gamma)*cos(chi)/(R*cos(lambda));
y(2) = V*cos(gamma)*sin(chi)/(R);
y(3) = V*sin(gamma);
y(4) = -g*sin(gamma) + Yfo/m + Omega^(2)*R*cos(lambda)*...
    (sin(gamma)*cos(lambda) - cos(gamma)*sin(chi)*sin(lambda));
y(5) = (V/R-g/V)*cos(gamma) + Xfo/(V*m) + 2*Omega*cos(chi)*cos(lambda)+ ...
    (Omega^2*R*cos(lambda)*(cos(gamma)*cos(lambda)+sin(gamma)*sin(chi)*sin(lambda)))/V;
y(6) = -V*cos(gamma)*cos(chi)*tan(lambda)/R + Zfo/(m*V*cos(gamma)) + ...
    2*Omega*(tan(gamma)*sin(chi)*cos(lambda)-sin(lambda))-(Omega^2*w3*cos(chi)*cos(lambda)*sin(lambda))/(V*cos(gamma));

y = y';
