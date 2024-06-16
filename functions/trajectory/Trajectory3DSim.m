function [simParams, trajectory] = Trajectory3DSim(x, density_mode, simInputs)

BC = x(1);
gamma0 = x(2);
mu = simInputs.mu;

params = [BC, density_mode, simInputs.R, simInputs.SC.m, simInputs.mu, simInputs.SC.L_D, simInputs.Omega];
x0 = [simInputs.delta0, simInputs.lambda0, simInputs.R0 , simInputs.V0, gamma0, simInputs.chi0];
[t, y] = ode45( @(t, x) reentryModel(t, x, params, simInputs.atmoModel), simInputs.tspan, x0, simInputs.odeOptions);

delta = y(:,1);     % Longitud
lambda = y(:,2);    % Latitude
r = y(:,3);         % Radial position
v = y(:,4);         % Velocity
gamma = y(:,5);     % FPA
chi = y(:,6);       % Heading angle
% an = axial_acceleration([r, v, gamma, lambda, chi], BC, simInputs);     % Axial acceleration
% [q_dot, q_dot_max, Q, Tw] = ThermalModel(simInputs, BC, [t,r,v], density_mode);
% trajectory = [t, r, v, gamma, delta, lambda, chi, q_dot, Q, an];

[q_dot, q_dot_max, Q, Tw] = ThermalModel(simInputs, BC, [t,r,v], density_mode);
trajectory = [t, r, v, gamma, delta, lambda, chi, q_dot, Q];

[r_vect, v_vect] = getECIframe(r, delta, lambda, v, gamma, pi/2 - chi);
r = r_vect(end,:)';
v =v_vect(end,:)';

[a, e, i, Omega, omega, theta] = getCOE(r, v, mu); % Classical orbital elemenets (COE)
simParams = [a, e, i, q_dot_max, Q(end), Tw(end), BC, gamma0, Omega, omega, theta];