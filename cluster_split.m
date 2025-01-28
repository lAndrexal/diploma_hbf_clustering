function [cl_directions, cl_inds, cl_sizes] = cluster_split(channels, N_clusters, N_iter,cluster_type, ind2train, cl_d)

N_rx = size(channels,2);
N_scen = size(channels,1); 
N_used = size(channels,3);

% consistent_inds = 1:140;
consistent_inds = load("consistent_inds").consistent_inds;

ch_directions = load("sv_svd_RX1024.mat").singular_vecs(consistent_inds,:);
% ch_directions = load("sa_svd_RX1024.mat").singular_subarrays(consistent_inds,:);
cl_sizes = zeros(N_clusters,1);
ch_directions = normalize_rows(ch_directions);
if cl_d == 0 
    cl_directions = normalize_elements(randn(N_clusters,N_rx) + 1i*randn(N_clusters,N_rx));
else
    cl_directions = cl_d / norm(cl_d);
end
if cluster_type == 1 % k-means
    for iter = 1:max(1,N_iter)
        sim_table = zeros(N_scen,N_clusters);

        % Create similarity table
        for subarray_ind = 1:size(ind2train,2)
            sim_table = sim_table + abs(ch_directions(:,ind2train(:,subarray_ind)) * cl_directions(:,ind2train(:,subarray_ind))');
        end

        
%         sim_table = abs(ch_directions(:,ind2train(:,:)) * cl_directions(:,ind2train(:,:))');
        
        [~,cl_inds] = max(sim_table, [], 2);
        f_dist = @(X1, X2)complex_subarray_cos(X1,X2, ind2train);
      
        % Get cluster statistics, cluster center - eigvec of concatenation
        if N_iter == 0
            for cl = 1:N_clusters
                cl_sizes(cl) = sum(cl_inds==cl);
            end
        else
            for cl = 1:N_clusters
                cl_sizes(cl) = sum(cl_inds==cl);
                cluster_elements = channels(cl_inds == cl, :, :);
                for i = 1:64

                    [v,e_v] = eigs(cluster_elements(:,ind2train(:,i))'*cluster_elements(:,ind2train(:,i)),1);
                    
                    cl_directions(cl,ind2train(:,i)) = normalize_elements(v);
                end
%                 cluster_elements = squeeze(reshape(permute(cluster_elements,[1,3,2]), [cl_sizes(cl)*N_used, N_rx]));
%                 [v,~] = eigs(cluster_elements'*cluster_elements,1);
%                 cl_directions(cl,:) = v;
            end
        end
    end
else % hierarchical44

    % Custom distance for hybrid beamforming
    f_dist = @(X1, X2)complex_subarray_cos(X1,X2, ind2train);
    % f_dist = @(X1, X2)complex_subarray(X1,X2);
    % Next functions are standard
    distance_vec = pdist(ch_directions, f_dist);
    linkage_result = linkage(distance_vec, 'complete');
    cl_inds = cluster(linkage_result, 'MaxClust',N_clusters);
%     silhouette(ch_directions,cl_inds,f_dist)
    % Get channel statistics
    for cl = 1:N_clusters
        cl_sizes(cl) = sum(cl_inds==cl);
        cluster_elements = channels(cl_inds == cl, :, :);
        cluster_elements = squeeze(reshape(permute(cluster_elements,[1,3,2]), [cl_sizes(cl)*N_used, N_rx]));
        if 0
            [v,~] = eigs(cluster_elements'*cluster_elements,1);
            cl_directions(cl,:) = v;
        else
            for i = 1:64

                [v,~] = eigs(cluster_elements(:,ind2train(:,i))'*cluster_elements(:,ind2train(:,i)),1);
                
%                 cl_directions(cl,ind2train(:,i)) = v / sqrt(64);
                cl_directions(cl,ind2train(:,i)) = steering_from_eig(v) / sqrt(16);
            end
        end
        
%         cl_ch = ch_directions(cl_inds==cl,:);
%         cl_directions(cl,:) = cl_ch(1,:);
    end
end
end

function A = normalize_rows(A)
    for i = 1:size(A,1)
        A(i,:) = A(i,:) / norm(A(i,:));
    end
end

function A = normalize_elements(A)
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            A(i,j) = A(i,j) / norm(A(i,j));
        end
    end
end

function steering_matrix = steering_from_eig(E)
   
    Nvert = 4;
    Nhor = 4;
    Npol = 1;

    C=[];
    D=[];
    for sub_array_index=1:Npol

        A=squeeze(E(sub_array_index:end).');
        
        G=reshape(A,[],Nhor);
        
        % angle 1 correlations
        for i=1:Nvert
            C=[C conj(G(i,1:Nhor-1)).*G(i,2:Nhor)];
        end
        
        % angle 2 correlations
        for j=1:Nhor
            D=[D (conj(G(1:Nvert-1,j)).*G(2:Nvert,j)).'];
        end

    end
    
    % //////////test antennas correlation////////////
    % cor_vec=(abs(corr_sum1)+abs(corr_sum2))/a(sample_indx);
    
    alpha=angle(sum(C));
    beta=angle(sum(D));
    
    steering_matrix = zeros(Nvert,Nhor);
        
    % make steering vector (the same for each sub-array)
    for k=1:Nvert
        steering_matrix(k,:)=exp(1i*(alpha*(0:1:Nhor-1).'+beta*(k-1)));
    end  
    steering_matrix = steering_matrix(:);
end