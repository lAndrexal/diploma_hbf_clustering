rng(42)
fig_ind = 1;
N_scen = 140;
N_digital = 64;
N_rx = 1024;
N_sc = 192;
N_ant = N_rx / N_digital;


dir = "./";
consistent_inds = load(dir + "consistent_inds").consistent_inds; % good LOS
% consistent_inds = 1:140;
% H = load(dir +"data_full_RX1024_split.mat").H_full;
H = load(dir +"data_full_RX1024_split.mat").H_full(consistent_inds,:,:);
N_scen = size(H,1);

inds_1D = 1:N_rx;
inds_2D = reshape(inds_1D, [2 16 32]);
ind2train = zeros(N_ant,N_digital);

for i = 1:N_digital

    [p,v,h] = ind2sub([2 4 8], i);
    sv = 4;
    sh = 4;
    ind = inds_2D(p, (v-1)*sv+1:v*sv, (h-1)*sh+1:h*sh);
    ind2train(:, i) = ind(:);
end

%%
% phases{1} = load(dir+"Trained_phases_uniform").phases;
% phases{2} = load(dir+"trained_phases_GA2.mat").phases;
% phases{3} = rand(1,N_rx) * 2 * pi;
% phases{4} = load(dir+"Analog_learningGD_UNIFORM.mat").phases;
% phases{5} = load(dir+"Analog_learningGD_UNIFORM_2.mat").phases;
if ~all(size(H) == [N_scen, N_rx, N_sc])
    warning("size of H is non consistent with predefined values N_")
end
% H_reshaped = squeeze(reshape(permute(H(1:N_scen, 1:N_rx, 1:N_sc),[2,1,3]), [N_rx, N_scen*N_sc]));
% [V,D] = eigs(H_reshaped*H_reshaped',1);
% phases{6} = angle(V);

norm_factor = sqrt(N_ant);


approaches_to_analyze = [0 -2 -4];
loss_ref = 1;
N_cluster = 16;
loss = zeros(numel(approaches_to_analyze),N_cluster,N_scen);
legend_entity = [];
approach_ind = 0;
% U = load("sa_svd_RX1024.mat").singular_subarrays(consistent_inds,:);
U = load("sv_svd_RX1024.mat").singular_vecs(consistent_inds,:);
for j = approaches_to_analyze 
    approach_ind = approach_ind + 1;
    for n_cl = N_cluster
        
    
       
        legend_entity = [legend_entity, "mode "+num2str(j)+", n_cl = "+num2str(n_cl)];
        
        if ismember(j, [1:6])
            for i = 1:N_scen
                phaseshifts = exp(1i*phases{j}) / norm_factor;
                loss(approach_ind,n_cl,i) = phaseshift_loss(phaseshifts,squeeze(H(i,:,:)), ind2train, 0);
            end
        elseif j == 10
            for i = 1:N_scen
                u = U(i,:).';
                
                % Here Fanalog has only 1 nonzero column (no blocks of subarrays)
                Fanalog = zeros(N_rx, N_digital);
                Fanalog(:, 1) = u;
                P = Fanalog'*squeeze(H(i,:,:));
                loss(approach_ind,n_cl,i) = - sum(abs(P).^2,'all');
    
            end
        elseif ismember(j,[-2 -1])
            
            [cl_d, cl_i, cl_s] = cluster_split(H,n_cl,20,abs(j), ind2train, 0);
            
            for i = 1:N_scen
                phaseshifts = exp(1i* angle(cl_d(cl_i(i),:)'))/norm_factor;
                loss(approach_ind,n_cl,i) = phaseshift_loss(phaseshifts,squeeze(H(i,:,:)), ind2train, 0);
            end

        elseif j == -4
          
            [V,D] = eigs(H_reshaped*H_reshaped',16);
            cl_d = V';

            [cl_d, cl_i, cl_s] = cluster_split(H,16,0,1, ind2train,cl_d);

            for i = 1:N_scen
                phaseshifts = exp(1i* angle(cl_d(cl_i(i),:)'))/norm_factor;
                loss(approach_ind,n_cl,i) = phaseshift_loss(phaseshifts,squeeze(H(i,:,:)), ind2train, 0);
            end
            
        elseif j == -3
            [cl_d, cl_i, cl_s] = cluster_split_R(H,n_cl,20,2, ind2train, 0);
            
            for i = 1:N_scen
                phaseshifts = exp(1i* angle(cl_d(cl_i(i),:)') ) /norm_factor;
                loss(approach_ind,n_cl,i) = phaseshift_loss(phaseshifts,squeeze(H(i,:,:)), ind2train, 0);
            end
            
        elseif j == 0
            for i = 1:N_scen
                u = U(i,:).';
                phaseshifts = exp(1i* angle(u) ) /norm_factor;
                % phaseshifts = u;
                loss(approach_ind,n_cl,i) = phaseshift_loss(phaseshifts,squeeze(H(i,:,:)), ind2train, 0);
            end
            loss_ref = loss(1,n_cl,:);
        end
    
        
        
        disp([j, n_cl, -sum(loss(approach_ind,n_cl,:),3)])
        if approach_ind == 1
            break;
        end
    end
end

% if approach_ind ~=1
figure(fig_ind)
hold on

%%
loss_ref = squeeze(loss_ref(1,1,:));
loss = squeeze(loss(:,N_cluster,:)).';
if 1
    % scatter(1:N_scen, loss ./ loss_ref);
    plot(1:N_scen, loss);
else
    boxplot(loss ./ loss_ref);
end
legend(legend_entity)