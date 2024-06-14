function [AllGamma, Params, Trajectories, DesignDrivers] = aeroCorridor(simInputs)

BC_range = simInputs.BC_range;
gamma_range = simInputs.gamma_range;
robustCorridor = simInputs.robustCorridor;

nBC = length(BC_range);

%% INITIALISING PARAMETERS
density_label = {'normal', 'low', 'high'};
if robustCorridor
    ndensity = [1,2,3];
    Params = zeros([nBC, ndensity*9*2]);
    DesignDrivers = zeros([nBC, ndensity*7*2]);
    AllGamma = zeros([nBC, 6]);
    Trajectories = cell([nBC, 1]);
end

switch simInputs.AerocaptureCorridor
    case 1  % ONLY RE-ENTRY TRAJECTORY
        Params = zeros([nBC, 9]);
        DesignDrivers = zeros([nBC, 7]);
        AllGamma = zeros([nBC, 6]);
        Trajectories = cell([nBC, 1]);
        for ibc = 1:nBC % Looping BC range
            fprintf("\n\tRunning Trajectory case:  (BC, FPA) = (%6.2f, %2.3f deg)", BC_range(ibc), gamma_range(ibc)*180/pi)

            [simParams, trajResults] = Trajectory3DSim([BC_range(ibc),gamma_range(ibc)], simInputs.densityMode, simInputs);
            t = trajResults(:,1);
            tend =  t(end);
            tq = linspace(0, tend, simInputs.nPoints);
            V = trajResults(:,2:end);
            Vq = interp1(t, V, tq, 'pchip');
            % F = griddedInterpolant(trajResults(:,1), V);
            % Vq = F(tq);
            array = [tq', Vq];

            [~, drivers] = cost_function(simInputs, simParams, simInputs.Opti);

            % OUTPUT
            Params(ibc,:) = [simParams(1), simParams(2), simParams(3), simParams(9), simParams(10), simParams(11), simParams(4), simParams(5), simParams(6)];
            DesignDrivers(ibc,:) = drivers;
            Trajectories{ibc} = array;
            AllGamma(ibc,:) = simParams(9);
        end
    case 2  % AEROCAPTURE CORRIDOR SERIAL IMPLEMENTATION
        if robustCorridor
            [gamma_lower, gamma_upper] = getRobustCorridor(simInputs, BC);
            Dgamma= (gamma_upper - gamma_lower)/2;

            [simParams, trajResults] = Trajectory3DSim([BC,gamma_lower], density_mode, simInputs);
            t = trajResults(:,1);
            tq = linspace(0, t(end), simInputs.nPoints);
            V = trajResults(:,2:end);
            F = griddedInterpolant(t, V);
            Vq = F(tq);
        else
            denMod = simInputs.densityMode;
        end
        for ibc = 1:nBC % Looping BC range
            tempGamma = zeros([1,2*3]);
            tempParams = zeros([1,2*3*9]);
            tempTraj = zeros([simInputs.nPoints, 3*2*9]);
            tempDrivers = zeros([1,2*3*7]);
            BC = BC_range(ibc);
            fprintf("\nSTARTING BC case " + ibc + " of " + nBC + " (" + BC + ").");

            for denMod = 1:N
                density_mode = ndensity(denMod);
                results = zeros([9,2]);
                drivers = zeros([7,2]);
                fprintf("\n\tDensity scenario " + density_mode + ": " + density_label{density_mode} + ".")
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
                [~, temp] = cost_function(simInputs, simParams, simInputs.Opti);
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
                [~, temp] = cost_function(simInputs, simParams, simInputs.Opti);
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
    case 3  % AEROCAPTURE CORRIDOR
        %% START PARALLEL COMPUTING POOL
        parfor ibc = 1:nBC % Looping BC range
            localSimInputs = simInputs;

            tempGamma = zeros([1,2*3]);
            tempParams = zeros([1,2*3*9]);
            tempTraj = zeros([localSimInputs.nPoints, 3*2*10]);
            tempDrivers = zeros([1,2*3*6])
            BC = BC_range(ibc);

            for density_mode = 1:ndensity
                results = zeros([9,2]);
                drivers = zeros([7,2])

                % Get corridor boundaries
                [gamma_lower, gamma_upper] = getCorridorBoundary(simInputs, BC, density_mode)
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
                %                 [~, temp] = cost_function(simInputs, simParams, simInputs.Opti);
                [~, temp] = cost_function(simInputs, simParams, simInputs.Opti);
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
                [~, temp] = cost_function(simInputs, simParams, simInputs.Opti);
                Dgamma = - Dgamma;
                drivers(:,2) = [temp(1:end-1),Dgamma];

                tempGamma(1,2*(density_mode-1)+1:2*density_mode) = [gamma_lower, gamma_upper];
                tempParams(2*9*(density_mode-1)+1:2*9*density_mode) = results(:)';
                tempTraj(:,20*(density_mode-1)+1:20*density_mode) = [array1, array2];
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
fprintf('\n')