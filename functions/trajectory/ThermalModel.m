function [q_dot, q_dot_max, Q, Tw] = ThermalModel(simInputs, BC, trajectoryData, density_mode)

if nargin < 3
    density_mode = 1;
end

RB = sqrt(simInputs.SC.m / ( BC * simInputs.SC.CD * pi));   % Aershoell base radius
RN = simInputs.SC.RN_RB*RB;

time = trajectoryData(:,1);
h = trajectoryData(:,2) - simInputs.R;
V = trajectoryData(:,3);
[~, ~, rho] = simInputs.atmoModel(h, density_mode, 0);

% radiative heat transfer
q_rad = simInputs.constTS{1} .* RN.^simInputs.constTS{2}(V, rho) .* rho.^simInputs.constTS{3} .* simInputs.f_V(V); % radiative heat flux from fluid to wall (W/m^2)


switch simInputs.thermalModel
    case 1 % Sutton-Graves
        q_conv = simInputs.k .* sqrt(rho./RN) .* V.^3; % convective heat flux from fluid to wall (W/m^2)
    case 2 % Fay-Ridell
        error('Fay-Ridell Stagnation point heating has not yet been implemented. Please use Sutton-Graves.')
    otherwise
        error('Please input a valid stagnation point thermal model index.')
end

q_dot = q_conv + q_rad; % net heat flux/unit area from fluid to wall (W/m^2)
Q = zeros(size(time));
for i = 2:length(time)
    Q(i) = trapz(time(1:i), q_dot(1:i)); % Heat load
end
q_dot_max = max(q_dot);
Tw = (q_dot_max./(simInputs.SC.epsilon*simInputs.sigma)).^(1/4);
