rng(42);

Nrx = 1024;
N_user = 4; % Number of users in simulation (8 - max)
N_used = 192;
N_seed = 8;

N_iter = 200; % Number of random selection of users to group


                    'Cluster correlation set';
% Filter out some users
consistent_inds = load("consistent_inds.mat").consistent_inds;     
% consistent_inds = 1:140;
N_scen = numel(consistent_inds);
precalculated_indices = zeros(N_scen,N_seed,N_user);

corr_mat = load("corrmat_svd_16RB_RX"+num2str(Nrx)).corr_mat(consistent_inds, consistent_inds); % Precalculated correlation matrix

cl_i = load("SVD16_check.mat").cl_i; % Load clusterization results (cl_d is used in beamspace_selection,m)

assert(all(corr_mat - diag(ones(N_scen,1)) < 1,'all')); % Check for identical users

th_max = 0.5;
for scen=1:N_scen % Iterate over all users (try to match everyone to some group)
    
    good_inds = find(cl_i == cl_i(scen))'; % All users from current cluster
    N_good = numel(good_inds);
    
    if N_good >= N_user 
        
        for seed = 1:N_seed % Try to create at least N_seed group with current user
            
            for i=1:N_iter % Randomly choose users from current cluster unless correlation threshold is met
                tmp_indexes=randperm(N_good-1,N_user-1);
                index = [scen good_inds(tmp_indexes)];
                sel_corr = corr_mat(index, index);

                max_corr = max(sel_corr -diag(ones(N_user,1)), [],'all'); % Delete 1s at diagonal
                if max_corr < th_max
                    precalculated_indices(scen, seed, :) = index;
                    break
                end
            end
        end
    end
end

% If some users were filtered out, match chosen ones with their original indices
for i = 1:N_scen
    for j = 1:N_seed
        if squeeze(precalculated_indices(i,j,:)) ~= zeros(N_user,1)
            precalculated_indices(i,j,:) = consistent_inds(squeeze(precalculated_indices(i,j,:)));
            
        end
    end
end