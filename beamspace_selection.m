function [beam_matrix]  = beamspace_selection(h_SRS, scenario)

% Initialize system parameters
N_ports=scenario.N_ports;
Ndigital = scenario.Ndigital;
Nrx = scenario.Nrx;
N_target_UE = scenario.N_target_UE;
training = scenario.training;
training_mode = scenario.training_mode;

% Initialize antenna lattice structure
N_pol_total = 2;
if scenario.Nrx == 64
    N_hor_total  = 8;
    N_vert_total = 4;
elseif scenario.Nrx == 128
    N_hor_total  = 8;
    N_vert_total = 8;
elseif scenario.Nrx == 256
    N_hor_total  = 16;
    N_vert_total = 8;
elseif scenario.Nrx == 1024
    N_hor_total  = 32;
    N_vert_total = 16;
end


if ~scenario.hybrid
    % Determine the antenna lattice parameters
    if scenario.beamspace_mode == 2
        Wfin  = (dftmtx(Nrx));
        D = abs(diag(Wfin'*Wfin));
        Wfin = Wfin * diag(1.0./sqrt(D));
    elseif scenario.beamspace_mode == -1
        if numel(size(h_SRS)) > 2
            h_SRS = reshape( permute( h_SRS ,[2 1 3]), Nrx, []);
        end
        [U,S,~] = svd(h_SRS);
        assert(isequal(sort(diag(S), "descend"),diag(S)))
        beam_matrix = U(:,1:N_ports);
        return
    end
    Nspace = Nrx;
else
    % h_SRS [Nue, Nrx, Nsc] 
    assert(scenario.Nrx == 1024);
    Wanalog = zeros(scenario.Nrx, scenario.Ndigital);
    
    % Create mapping between subarray and antennas
    inds_1D = 1:1024;
    inds_2D = reshape(inds_1D, [N_pol_total N_vert_total N_hor_total]);
    if scenario.digital_pol % Here can be other structures (e.g. 16-2-2)
        Nsa_vert = 4;
        Nsa_hor = 8;
        Nsa_pol = 2;
    else
        Nsa_vert = 8;
        Nsa_hor = 8;
        Nsa_pol = 1;
    end
    sa2rx_inds = zeros(Ndigital, Nrx/Ndigital);
    for i = 1:scenario.Ndigital
        [p,v,h] = ind2sub([Nsa_pol Nsa_vert Nsa_hor], i);
        N_hor_perSA = N_hor_total / Nsa_hor;
        N_vert_perSA = N_vert_total / Nsa_vert;
        ind = inds_2D(p, (v-1)*N_vert_perSA+1:v*N_vert_perSA, (h-1)*N_hor_perSA+1:h*N_hor_perSA);
        sa2rx_inds(i,:) = ind(:);
    end
    
    % Determine values at phaseshifters
    phases = zeros(1024,1);
    switch scenario.analog_mode
        case 0 % Just looking forward
            % same zeros
        case 1 % Offline trained on 140 scenarios
            phases = load("Beamforming/phases_analog.mat").phases;
        case 2 % Online trained for each scenario
            scen_index = (scenario.index-1)*16+scenario.seed; 
%             phases = gd_analog_call(h_SRS, scenario.Ndigital, load("Beamforming/Python_input.mat").ind2train);
%             save("Beamforming/saved_phases/phases_"+num2str(scen_index)+"scen", 'phases');
            phases = load("Beamforming/saved_phases/phases_"+num2str(scen_index)+"scen").phases;
        case 3 % Random phases (no channel info?)
            phases = (2*rand(1024,1)-1) * pi;
        case 4 % Trained with genetic algorithm
            phases = load("Beamforming/Trained_phases_GA2.mat").phases;
        case 5 %  Trained with uniformity term
            phases = load("Beamforming/Trained_phases_uniform.mat").phases;
        case 6 % Long trained with gradient descent and uniformity 
            phases = load("Beamforming/Analog_learningGD_UNIFORM_2.mat").phases;
        case -2 % For subarrays SU whole SVD phases
            [u,~,~]=svds(squeeze(h_SRS(1,:,:)),1);
            phases = angle(u);
        case -3 % For subarray set MU whole SVD phases
            Nsa_perUE = Ndigital / N_target_UE;
            for user_id = 0:N_target_UE-1
                [u,~,~]=svd(squeeze(h_SRS(user_id+1,:,:)));
                for i = int32(user_id*Nsa_perUE+1:(user_id+1)*Nsa_perUE)
                    phases(sa2rx_inds(i,:)) = angle(u(sa2rx_inds(i,:),1));
                end
            end
        case -4 % For subarray set MU subarray SVD phases
            Nsa_perUE = Ndigital / N_target_UE;
            for user_id = 0:scenario.N_target_UE-1
                size_subarray = Nrx / Ndigital;
                arrays_ind = user_id*Nsa_perUE+1:(user_id+1)*Nsa_perUE;
                antennas_inds = sa2rx_inds(arrays_ind,:);
                channel_part = squeeze(h_SRS(user_id+1,antennas_inds,:)); 
                [u,~,~]=svds(channel_part,1);
                u = reshape(u, [size_subarray,Nsa_perUE]);
                for i = 1:numel(arrays_ind) 
                    phases(antennas_inds(i,:)) = angle(u(i,:));
                end
            end
        case -5 % For subarrays MU whole SVD phases
            H = reshape(permute(h_SRS, [2,1,3]), scenario.Nrx, []);
            [u,~,~] = svds(H,1);
            phases = angle(u);
        case -6 % For subarrays MU cluster SVD phases
            u = load("Clustering&Learning/pi_subarray_svd_RX1024.mat").groups_directions(scenario.index, scenario.seed,:);
            phases = angle(squeeze(u));
        case -7 % For subarrays MU cluster SVD phases (filtered, GA)
            load('Clustering/consistent_inds', 'consistent_inds'); % channels with dominating LOS component

            % mapping between all and consistent channels and cluster index
            cl_i = load(scenario.clustering_path).cl_i(ismember(consistent_inds, scenario.precalculated_indices(scenario.index,scenario.seed,1)));
            
            % choosing the group direction
            u = load(scenario.clustering_path).cl_d(cl_i,:);
            if scenario.clustering_path == "Dergachev\Clustering&Learning\HC_sub"
                phases = -angle(squeeze(u));
            else
                phases = angle(squeeze(u));
            end

        case -8 % For subarrays MU whole steering phases
            steering = reshape(get_steering_matrix(squeeze(h_SRS)), [512,1]);
            steering_vector = zeros(1024,1);
            steering_vector(1:2:1024) = steering;
            steering_vector(2:2:1024) = steering;
            phases = angle(steering_vector);
        case -9 % For subarrays MU subarray steering phases     
            h_SRS2 = squeeze(h_SRS);
             for i = 1:scenario.Ndigital
                phases(sa2rx_inds(i,:) ) = angle( reshape( get_steering_matrix( h_SRS2(sa2rx_inds(i,:) )), [16,1]) );
             end
        case -10 % For subarrays SU subarray SVD values
            h_SRS2 = squeeze(h_SRS);
             for i = 1:scenario.Ndigital
                Wanalog(sa2rx_inds(i,:), i) = svds(h_SRS2(sa2rx_inds(i,:)), 1);
             end

        case -12 % For whole SU whole SVD values
            h_SRS2 = squeeze(h_SRS);
            [u, ~, ~] = svds(h_SRS2,1);
            Wanalog(:,1) = u;   
    end
    
    % Put computed phases to analog beamforming matrix
    if ~ismember(scenario.analog_mode, [-10, -12])
        for i = 1:scenario.Ndigital
            ind = sa2rx_inds(i,:);
            if scenario.subarray_phase
                norm_multiplier = 1;
            else
                norm_multiplier = exp(1i*phases(ind(1)));
            end
            Wanalog(ind(:), i) = exp(1i*phases(ind(:))) / norm_multiplier;
        end
    end

    % To have the right scale in detector
    if ismember(scenario.analog_mode, [1 -10:-1])
        norm_coefficient = sqrt(16);
    elseif scenario.analog_mode == -12
        norm_coefficient = 1/8;
    end
    Wanalog = Wanalog / norm_coefficient; 
    
    % Digital part beamforming
    switch scenario.beamspace_mode 
        case 2 % Simple DFT beams
            Wdigital  = (dftmtx(scenario.Ndigital));
            D = abs(diag(Wdigital'*Wdigital));
            Wdigital = Wdigital * diag(1.0./sqrt(D));
        case 4 % Offline trained FFT (with offline trained phases)
            Wdigital = load("WeightsHybrid_mode3.mat").T_N;
        case 12
            if training; [~, T_N] = gd_digital(h_SRS, Wanalog, training_mode, 32); end
            scen_index = (scenario.index-1)*16+scenario.seed; 
            if scenario.analog_mode == 4
                if training; save("Beamforming/saved_digital/offline_phases_digital_"+num2str(scen_index)+"scen", 'T_N'); end
                Wdigital = load("Beamforming/saved_digital/offline_phases_digital_"+num2str(scen_index)+"scen").T_N';
            elseif scenario.analog_mode == 2
                if training; save("Beamforming/saved_digital/online_phases_digital_"+num2str(scen_index)+"scen", 'T_N'); end
                Wdigital = load("Beamforming/saved_digital/online_phases_digital_"+num2str(scen_index)+"scen").T_N';
            elseif scenario.analog_mode == 3
                if training; save("Beamforming/saved_digital/random_phases_digital_"+num2str(scen_index)+"scen", 'T_N');end
                Wdigital = load("Beamforming/saved_digital/random_phases_digital_"+num2str(scen_index)+"scen").T_N';
            elseif scenario.analog_mode == 5
                if training; save("Beamforming/saved_digital/uniform_phases_digital_"+num2str(scen_index)+"scen_mode"+num2str(scenario.training_mode), 'T_N'); end
                Wdigital = load("Beamforming/saved_digital/uniform_phases_digital_"+num2str(scen_index)+"scen_mode"+num2str(scenario.training_mode)).T_N';
            elseif scenario.analog_mode == 6
                if training; save("Beamforming/saved_digital/uniformGD_phases_digital_"+num2str(scen_index)+"scen_mode"+num2str(scenario.training_mode), 'T_N');end
                Wdigital = load("Beamforming/saved_digital/uniformGD_phases_digital_"+num2str(scen_index)+"scen_mode"+num2str(scenario.training_mode)).T_N';
            elseif scenario.analog_mode == -7
                if training; save("Beamforming/saved_digital/filtered_GAcl_phases_digital_"+num2str(scen_index)+"scen_mode"+num2str(scenario.training_mode), 'T_N');end
                Wdigital = load("Beamforming/saved_digital/filtered_GAcl_phases_digital_"+num2str(scen_index)+"scen_mode"+num2str(scenario.training_mode)).T_N';
            end
        case -1
            if numel(size(h_SRS)) > 2
                h_SRS = reshape( permute( h_SRS ,[2 1 3]), scenario.Nrx, []);
            end
            h_SRS = Wanalog' * h_SRS;
            [U,S,~] = svd(h_SRS);
            assert(isequal(sort(diag(S), "descend"),diag(S)))
            beam_matrix = Wanalog * U(:,1:N_ports);
            return
        case 0
            Wdigital = eye(scenario.Ndigital);
    end
    Wfin = Wanalog*Wdigital;
    Nspace = Ndigital;
end

% Choosing most powerful beams among digital ones
beam_power_norm = zeros(scenario.N_user, Nspace);
for i = 1 : scenario.N_user
    if i <= N_target_UE
        beam_power = sum( abs( squeeze(h_SRS(i,:,:)).' * conj(Wfin)).^2 , 1);
    else
        beam_power = 1;
    end
    beam_power_norm(i,:) = beam_power / sum(beam_power);
end
beam_power_MU = sum (beam_power_norm,1);

[~, beam_indx] = sort(beam_power_MU, 'descend');
beam_matrix=Wfin(:,beam_indx(1:N_ports));