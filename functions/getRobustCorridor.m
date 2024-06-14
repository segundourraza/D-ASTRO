function [gamma_lower, gamma_upper] = getRobustCorridor(simInputs, BC)

Rplanet = simInputs.R;
h_AI = simInputs.h_AI;
mu = simInputs.mu;

% START BINRAY SEARCH FOR LOWER ROBUST AEROCAPTURE CORRIDOR BOUNDARY
% Corresponds to HIGH density mode
params = [BC, 3, simInputs.R, simInputs.SC.m, simInputs.mu, simInputs.SC.L_D, simInputs.Omega];

gamma_min = simInputs.gamma_min;
gamma_max = simInputs.gamma_max;
gamma_lower = 3*simInputs.tol;
gamma_prev = 0;
iteration = 1;

while abs(gamma_lower-gamma_prev) > simInputs.tol && (iteration < simInputs.MaxIterations)
    gamma_prev = gamma_lower;
    gamma_lower = 0.5*(gamma_max + gamma_min);
    x0 = [simInputs.delta0, simInputs.lambda0, simInputs.R0, simInputs.V0, gamma_lower, simInputs.chi0];
    [~, y] = ode45( @(t, x) reentryModel(t, x, params, simInputs.models), simInputs.tspan, x0, simInputs.Opt);
    rEnd = y(end,3); % Radius of SC at end of trajectory
    
    % Identify lower limit of corridor, separates surface collision with no surface collision
    if (rEnd - Rplanet) < h_AI
        gamma_min = gamma_lower;
    else
        gamma_max = gamma_lower;
    end
    iteration = iteration + 1;
end

% START BINRAY SEARCH FOR UPPER ROBUST AEROCAPTURE CORRIDOR BOUNDARY
% Corresponds to LOW density mode
params(2) = 2;
gamma_min = gamma_lower;
gamma_max = simInputs.gamma_max;
gamma_upper = 3*simInputs.tol;
gamma_prev = 0;
iteration = 1;

while abs(gamma_upper-gamma_prev) > simInputs.tol && (iteration < simInputs.MaxIterations)
    gamma_prev = gamma_upper;
    gamma_upper = 0.5*(gamma_max + gamma_min);
    x0(5) =  gamma_upper;
    [~, y] = ode45( @(t, x) reentryModel(t, x, params, simInputs.models), simInputs.tspan, x0, simInputs.Opt);
    rEnd = y(end,3); % Altitude of SC at end of trajectory
    vEnd = y(end,4); % Velocity of SC at end of trajectory
    specific_energy = 0.5*vEnd^2 - mu/rEnd;
    % Identify lower limit of corridor, separated by negative orbital 
    % specific
    if specific_energy < 0
        gamma_min = gamma_upper;
    else
        gamma_max = gamma_upper;
    end
    iteration = iteration + 1;
end