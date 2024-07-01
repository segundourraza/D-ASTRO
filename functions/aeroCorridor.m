function [AllGamma, Params, Trajectories, DesignDrivers] = aeroCorridor(simInputs)

BC_range = simInputs.BC_range;
nBC = length(BC_range);

%% INITIALISING PARAMETERS
density_label = {'normal', 'low', 'high'};

switch simInputs.AerocaptureCorridor
    case 1  % ONLY RE-ENTRY TRAJECTORY
        Params = zeros([nBC, 9]);
        DesignDrivers = zeros([nBC, 7]);
        AllGamma = zeros([nBC, 6]);
        Trajectories = cell([nBC, 1]);
        for ibc = 1:nBC % Looping BC range
            fprintf("\n\tRunning Trajectory case:  (BC, FPA) = (%6.2f, %2.3f deg)", BC_range(ibc), simInputs.gamma_range(ibc)*180/pi)

            [simParams, trajResults] = Trajectory3DSim([BC_range(ibc),simInputs.gamma_range(ibc)], simInputs.densityMode, simInputs);
            t = trajResults(:,1);
            tend =  t(end);
            tq = linspace(0, tend, simInputs.nPoints);
            Vq = interp1(t, trajResults(:,2:end), tq, 'pchip');

            % OUTPUT
            Params(ibc,:) = [simParams(1), simParams(2), simParams(3), simParams(9), simParams(10), simParams(11), simParams(4), simParams(5), simParams(6)];
            [~, DesignDrivers(ibc,:)] = cost_function(simInputs, simParams);
            Trajectories{ibc} = [tq', Vq];
            AllGamma(ibc,:) = simParams(9);
        end
    case 2  % AEROCAPTURE CORRIDOR SERIAL IMPLEMENTATION
        if  ~simInputs.completeCorridor
            Params = zeros([nBC, 9*2]);             % Stores trajectory parameters
            DesignDrivers = zeros([nBC, 7*2]);      % Stores design driver values
            AllGamma = zeros([nBC, 2]);             % Stores gammas at boundaries
            Trajectories = [];                      % Stores trajctories
            for ibc = 1:nBC
                fprintf("\n\tStrating BC case %4u of %4u (%6.2f)", ibc, nBC, BC_range(ibc));
                % Compute lower and upper robust FPA angles (binary search)
                [AllGamma(ibc,1), AllGamma(ibc,2)] = getRobustCorridor(simInputs, BC_range(ibc));

                Dgamma= (AllGamma(ibc,2) - AllGamma(ibc,1))/2;

                % Running LOWER boundary trajectory with HIGH density mode
                % to get driver values
                [simParams, ~] = Trajectory3DSim([BC_range(ibc),AllGamma(ibc,1)], 3, simInputs);
                [~, temp] = cost_function(simInputs, simParams);
                % Add simualtion results:
                % [a, e, i, Omega, omega, theta,qdot_max, Q, Tw]
                Params(ibc,1:9) = [simParams(1), simParams(2), simParams(3), simParams(9), simParams(10), simParams(11), simParams(4), simParams(5), simParams(6)];
                DesignDrivers(ibc,1:7) = [temp(1:end-1), Dgamma];

                % Running UPPER boundary trajectory with LOW density mode
                % to get driver values
                [simParams, ~] = Trajectory3DSim([BC_range(ibc),AllGamma(ibc,2)], 2, simInputs);
                [~, temp] = cost_function(simInputs, simParams);
                % Add simualtion results:
                % [a, e, i, Omega, omega, theta,qdot_max, Q, Tw]
                Params(ibc,10:18) = [simParams(1), simParams(2), simParams(3), simParams(9), simParams(10), simParams(11), simParams(4), simParams(5), simParams(6)];
                DesignDrivers(ibc,8:14) = [temp(1:end-1), Dgamma];
            end
        else
            Params = zeros([nBC, 3*9*2]);             % Stores trajectory parameters
            DesignDrivers = zeros([nBC, 3*7*2]);      % Stores design driver values
            AllGamma = zeros([nBC, 3*2]);             % Stores gammas at boundaries
            Trajectories = cell(nBC);                 % Stores trajctories
            for ibc = 1:nBC % Looping BC range
                tempGamma = zeros([1,2*3]);
                tempParams = zeros([1,2*3*9]);
                tempTraj = zeros([simInputs.nPoints, 3*2*9]);
                tempDrivers = zeros([1,2*3*7]);
                BC = BC_range(ibc);
                fprintf("\n\tStrating BC case %4u of %4u (%6.2f)", ibc, nBC, BC_range(ibc));

                for density_mode = 1:3
                    results = zeros([9,2]);
                    drivers = zeros([7,2]);
                    fprintf("\n\t\tDensity scenario " + density_mode + ": " + density_label{density_mode} + ".")
                    % Get corridor boundaries
                    [gamma_lower, gamma_upper] = getCorridorBoundary(simInputs, BC, density_mode);
                    Dgamma= (gamma_upper - gamma_lower)/2;

                    [simParams, trajResults] = Trajectory3DSim([BC,gamma_lower], density_mode, simInputs);

                    t = trajResults(:,1);
                    tq = linspace(0, t(end), simInputs.nPoints);
                    V = trajResults(:,2:end);
                    F = griddedInterpolant(t, V);
                    Vq = F(tq);

                    array1 = [tq', Vq];
                    results(1,1) = simParams(1); % a
                    results(2,1) = simParams(2); % e
                    results(3,1) = simParams(3); % i
                    results(4,1) = simParams(9); % Omega
                    results(5,1) = simParams(10); % omega
                    results(6,1) = simParams(11); % theta
                    results(7,1) = simParams(4); % q_dot_max
                    results(8,1) = simParams(5); % Q(end)
                    results(9,1) = simParams(6); % Tw(end)
                    [~, temp] = cost_function(simInputs, simParams);
                    drivers(:,1) = [temp(1:end-1),Dgamma];

                    [simParams, trajResults] = Trajectory3DSim([BC,gamma_upper], density_mode, simInputs);

                    t = trajResults(:,1);
                    tq = linspace(0, t(end), simInputs.nPoints);
                    V = trajResults(:,2:end);
                    F = griddedInterpolant(t, V);
                    Vq = F(tq);

                    array2 = [tq', Vq];
                    results(1,2) = simParams(1);
                    results(2,2) = simParams(2);
                    results(3,2) = simParams(3);
                    results(4,2) = simParams(9);
                    results(5,2) = simParams(10);
                    results(6,2) = simParams(11);
                    results(7,2) = simParams(4);
                    results(8,2) = simParams(5);
                    results(9,2) = simParams(6);
                    [~, temp] = cost_function(simInputs, simParams);
                    Dgamma = (gamma_lower- gamma_upper)/2;
                    drivers(:,2) = [temp(1:end-1),Dgamma];

                    tempGamma(1,2*(density_mode-1)+1:2*density_mode) = [gamma_lower, gamma_upper];
                    tempParams(2*9*(density_mode-1)+1:2*9*density_mode) = results(:)';
                    tempTraj(:,18*(density_mode-1)+1:18*density_mode) = [array1, array2];
                    tempDrivers(:,2*7*(density_mode-1)+1:2*7*density_mode) = drivers(:)';
                end
                %% OUTPUT
                Params(ibc, :) = tempParams;
                Trajectories{ibc, 1} = tempTraj;
                AllGamma(ibc, :) = tempGamma;
                DesignDrivers(ibc,:) = tempDrivers;
            end
        end
    case 3  % AEROCAPTURE CORRIDOR
        %% START PARALLEL COMPUTING POOL

        if  ~simInputs.completeCorridor
            Params = zeros([nBC, 9*2]);             % Stores trajectory parameters
            DesignDrivers = zeros([nBC, 7*2]);      % Stores design driver values
            AllGamma = zeros([nBC, 2]);             % Stores gammas at boundaries
            Trajectories = [];                      % Stores trajctories
            fprintf("\n")
            parfor ibc = 1:nBC
                localSimInputs = simInputs;
                tempParams = zeros([1, 9*2]);
                tempDrivers = zeros([1, 7*2]);
                tempGamma = zeros([1, 2]);

                fprintf("\tStrating BC case %4u of %4u (%6.2f)\n", ibc, nBC, BC_range(ibc));
                % Compute lower and upper robust FPA angles (binary search)
                [tempGamma(1), tempGamma(2)] = getRobustCorridor(simInputs, BC_range(ibc));

                Dgamma= (tempGamma(2) - tempGamma(1))/2;

                % Running LOWER boundary trajectory with HIGH density mode
                % to get driver values
                [simParams, ~] = Trajectory3DSim([BC_range(ibc),tempGamma(1)], 3, simInputs);
                [~, temp] = cost_function(simInputs, simParams);
                % Add simualtion results:
                % [a, e, i, Omega, omega, theta,qdot_max, Q, Tw]
                tempParams(1:9) = [simParams(1), simParams(2), simParams(3), simParams(9), simParams(10), simParams(11), simParams(4), simParams(5), simParams(6)];
                tempDrivers(1:7) = [temp(1:end-1), Dgamma];

                % Running UPPER boundary trajectory with LOW density mode
                % to get driver values
                [simParams, ~] = Trajectory3DSim([BC_range(ibc),tempGamma(2)], 2, simInputs);
                [~, temp] = cost_function(simInputs, simParams);
                % Add simualtion results:
                % [a, e, i, Omega, omega, theta,qdot_max, Q, Tw]
                tempParams(10:18) = [simParams(1), simParams(2), simParams(3), simParams(9), simParams(10), simParams(11), simParams(4), simParams(5), simParams(6)];
                tempDrivers(8:14) = [temp(1:end-1), Dgamma];

                Params(ibc,:) = tempParams;
                DesignDrivers(ibc,:) = tempDrivers;
                AllGamma(ibc,:) = tempGamma;
            end
        else
            Params = zeros([nBC, 3*9*2]);             % Stores trajectory parameters
            DesignDrivers = zeros([nBC, 3*7*2]);      % Stores design driver values
            AllGamma = zeros([nBC, 3*2]);             % Stores gammas at boundaries
            Trajectories = cell([nBC, 2]);          % Stores trajctories
            parfor ibc = 1:nBC % Looping BC range
                localSimInputs = simInputs;

                tempGamma = zeros([1,2*3]);
                tempParams = zeros([1,2*3*9]);
                tempTraj = zeros([localSimInputs.nPoints, 3*2*10]);
                tempDrivers = zeros([1,2*3*6]);
                BC = BC_range(ibc);

                for density_mode = 1:3
                    results = zeros([9,2]);
                    drivers = zeros([7,2]);

                    % Get corridor boundaries
                    [gamma_lower, gamma_upper] = getCorridorBoundary(simInputs, BC, density_mode);
                    Dgamma = (gamma_upper - gamma_lower)/2;

                    [simParams, trajResults] = Trajectory3DSim([BC,gamma_lower], density_mode, simInputs);

                    t = trajResults(:,1);
                    tq = linspace(0, t(end), simInputs.nPoints);
                    V = trajResults(:,2:end);
                    F = griddedInterpolant(t, V);
                    Vq = F(tq);

                    array1 = [tq', Vq];
                    results(1,1) = simParams(1);
                    results(2,1) = simParams(2);
                    results(3,1) = simParams(3);
                    results(4,1) = simParams(9);
                    results(5,1) = simParams(10);
                    results(6,1) = simParams(11);
                    results(7,1) = simParams(4);
                    results(8,1) = simParams(5);
                    results(9,1) = simParams(6);
                    [~, temp] = cost_function(simInputs, simParams);
                    drivers(:,1) = [temp(1:end-1),Dgamma];

                    [simParams, trajResults] = Trajectory3DSim([BC,gamma_upper], density_mode, simInputs);

                    t = trajResults(:,1);
                    tq = linspace(0, t(end), simInputs.nPoints);
                    V = trajResults(:,2:end);
                    F = griddedInterpolant(t, V);
                    Vq = F(tq);

                    array2 = [tq', Vq];
                    results(1,2) = simParams(1);
                    results(2,2) = simParams(2);
                    results(3,2) = simParams(3);
                    results(4,2) = simParams(9);
                    results(5,2) = simParams(10);
                    results(6,2) = simParams(11);
                    results(7,2) = simParams(4);
                    results(8,2) = simParams(5);
                    results(9,2) = simParams(6);
                    [~, temp] = cost_function(simInputs, simParams);
                    Dgamma = - Dgamma;
                    drivers(:,2) = [temp(1:end-1),Dgamma];

                    tempGamma(1,2*(density_mode-1)+1:2*density_mode) = [gamma_lower, gamma_upper];
                    tempParams(2*9*(density_mode-1)+1:2*9*density_mode) = results(:)';
                    tempTraj(:,18*(density_mode-1)+1:18*density_mode) = [array1, array2];
                    tempDrivers(:,2*7*(density_mode-1)+1:2*7*density_mode) = drivers(:)';
                end

                %% OUTPUT
                Params(ibc, :) = tempParams;
                Trajectories{ibc, 1} = tempTraj;
                AllGamma(ibc, :) = tempGamma;
                DesignDrivers(ibc,:) = tempDrivers;

                fprintf("\tBC iter: " + ibc + "/" + nBC + " DONE.\n");
            end
        end
end
fprintf('\n')