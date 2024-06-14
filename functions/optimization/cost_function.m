function [cost, x] = cost_function(simInputs, simParams, Normalise)
if nargin < 3
    Normalise = 0;
end

%% CONSTANTS

% Planetary constants
mu = simInputs.mu;

% Trajectory Results
a = simParams(1); % Semi-major axis of post aerocapture orbit
e =simParams(2); % Eccentricity of post aerocapture orbit
i = simParams(3); % Inclination of post aerocapture orbit
q_dot_max = simParams(4); % Maximum heat flux of aerocapture trajectory
Q = simParams(5); % Heat load of aerocapture trajectory
Tw = simParams(6); % Rough wall temperature
BC = simParams(7);
gamma = simParams(8);

i_target = simInputs.Opti.i_orb_target;
a_target = simInputs.Opti.a_orb_target;
e_target = simInputs.Opti.e_orb_target;
%% PERFORMANCE DRIVERS
% Dimensional
V_umbrella = sqrt(simInputs.SC.m / ( BC * simInputs.SC.CD * pi));

% DELTA V COST FUNCTION
DV_total = getCorectiveManeuvreBurn(a, e, i, a_target, i_target, mu, simInputs.R0);

% DISTANCE FROM CORRIDOR BISECTION
if isfield(simInputs, 'fitlb')
    Dgamma = min(abs([simInputs.fitlb(BC),simInputs.fitub(BC)]- gamma));
else
    Dgamma = 0;
end

x = [V_umbrella, DV_total, Q, q_dot_max, Tw, a ,Dgamma];
%% FEATURE SCALING AND COST
cost = 0;
if Normalise
    x_norm = (x-simInputs.Opti.mins)./(simInputs.Opti.maxs -simInputs.Opti.mins); % Normalisation
    if simInputs.Opti.costFunc == 1
        cost = least_square_method(x_norm, simInputs.Opti.target, simInputs.Opti.weights);
    elseif simInputs.Opti.costFunc == 2
        cost = pseudo_huber_loss(x_norm, simInputs.Opti.target, simInputs.Opti.weights, simInputs.Opti.delta) ;
    end
end

end
