function [h_pilot_out, h_data_out, h_srs_out, UE_index]  = get_channel(scenario)

N_user=length(scenario.UE_power); % number of target users in scenario
Nrx = scenario.Nrx;               % N-of-Rx antennas
RB_size=scenario.RB_size;         % RB_size = 12 subcarriers
RB_num=scenario.RB_num;    
N_used=RB_num*RB_size;            % N-of-subcarriers (max=600)  
N_data_sym=12;
N_pilot_sym = 2;                  % N-of-pilot-symbols
N_target=scenario.N_target_UE;    % number of Interf users in scenario

% N_ofdm = 14;                      % N of symbols per TTI

% Data positions in TTI
data_positions1=1:3;
data_positions2=5:10;
data_positions3=12:14;

pilot_positions=[4 11]; % Pilot positions in TTI

h_pilot = zeros(N_user,Nrx, N_used,2);
h_data = zeros(N_user,Nrx, N_used,N_data_sym);

% If precalculated indices are available, fast loading
if ~isfield(scenario, 'precalculated_indices') 
    
    N_target=scenario.N_target_UE;       % number of Interf users in scenario
    N_scenarios=scenario.N_scenarios;
   
    index=num2str(scenario.index);
    
    % Load channel for the first user with the FIXED index !!!
    H_tmp = load_channel(index, Nrx, 0);

    h_pilot=zeros(N_user, Nrx,N_used,N_pilot_sym);
    h_data=zeros(N_user, Nrx,N_used,N_data_sym);
    
    comb = scenario.comb;
    
    % Extract pilot symbols
    h_pilot(1,:,1:(1+comb):N_used,1)= (sqrt(1+comb))*H_tmp(1,1:Nrx,1:(1+comb):N_used, pilot_positions(1));
    h_pilot(1,:,1:(1+comb):N_used,2) = (sqrt(1+comb))*H_tmp(1,1:Nrx,1:(1+comb):N_used,pilot_positions(2));

  % Extract channel for data symbols
    h_data(1,:,1:N_used, data_positions1-0) = H_tmp(1,1:Nrx,1:N_used, data_positions1);
    h_data(1,:,1:N_used, data_positions2-1) = H_tmp(1,1:Nrx,1:N_used, data_positions2);
    h_data(1,:,1:N_used, data_positions3-2) = H_tmp(1,1:Nrx,1:N_used, data_positions3);
    
    UE_index = zeros(1,N_user);
    UE_index(1)=scenario.index;
    
    % if the number of users (target+interf) > 1 then generate more
    if N_user > 1
        
        MU_index_set=[1:scenario.index-1 scenario.index+1:N_scenarios];   % exclude 1-st UE scenario
    
        if scenario.corr_flag == 1
            corr_tmp = 0;
        end
        
        if scenario.corr_flag == 0
            corr_tmp = 1;
        end
        
        for i=1:scenario.N_iter
            
            tmp_indexes=randperm(N_scenarios-1,N_user-1); % generate random indexes
            
            for j=2 : N_user
                
                if j<=N_target
                    comb = scenario.comb;
                else
                    comb = scenario.int_comb;
                end
                
                UE_index(j)=MU_index_set(tmp_indexes(j-1));
                index=num2str(UE_index(j));                
                H_tmp = load_channel(index, Nrx, 0);
                
                first_ind = 1;
                if comb && scenario.comb_split
                    first_ind = 2-mod(j,2); % 1 for odd, 2 for even
                end
                
                % Extract channel for pilot symbols
                h_pilot(j,:,first_ind:(1+comb):N_used,1) = (sqrt(1+comb))*H_tmp(1,1:Nrx,first_ind:(1+comb):N_used, pilot_positions(1));
                h_pilot(j,:,first_ind:(1+comb):N_used,2) = (sqrt(1+comb))*H_tmp(1,1:Nrx,first_ind:(1+comb):N_used, pilot_positions(2));
                
%                 % Extract channel for pilot symbols
%                 h_pilot(j,:,1:(1+comb):N_used,1) = (sqrt(1+comb))*H_tmp(1,1:Nrx,1:(1+comb):N_used, pilot_positions(1));
%                 h_pilot(j,:,1:(1+comb):N_used,2) = (sqrt(1+comb))*H_tmp(1,1:Nrx,1:(1+comb):N_used, pilot_positions(2));
                
                % Extract channel for data symbols
                h_data(j,:,1:N_used, data_positions1-0) = H_tmp(1,1:Nrx,1:N_used, data_positions1);
                h_data(j,:,1:N_used, data_positions2-1) = H_tmp(1,1:Nrx,1:N_used, data_positions2);
                h_data(j,:,1:N_used, data_positions3-2) = H_tmp(1,1:Nrx,1:N_used, data_positions3);
    
            end
            
            RB_num_new=RB_num/scenario.RB_num_Ruu;
            RB_size_new =RB_size*scenario.RB_num_Ruu;
            % test correlation
            counter=0;
            corr = zeros(1,N_user*(N_user-1) /2);
            corr_RB = zeros(RB_num_new, Nrx, Nrx);
            for m=1:N_user                 %  all users
                for n=m+1:N_user           %  all users
                    counter=counter+1;    % number of UE pairs
                    for k=1:RB_num_new           % calculate correlation for each RB
                        A1=squeeze(h_data(m,:,(k-1)*RB_size_new+1:k*RB_size_new,1)).';
                        A2=squeeze(h_data(n,:,(k-1)*RB_size_new+1:k*RB_size_new,1)).';
                        B1=diag(A1.'*conj(A1));
                        B2=diag(A2.'*conj(A2));
                        corr_RB(k,:,:)=abs(A1.'*conj(A2)) ./  (mean(sqrt(B1).*sqrt(B2)));  % correlation per RB
                    end
                    corr(counter)=mean(corr_RB,'all');  % average correlation between users i and j
                end
            end
            %[i min(corr) max(corr)]

             if (scenario.corr_flag ~= 1) && (scenario.corr_flag ~= 0)
                break
             end
            
            
             if scenario.corr_flag == 1
                 if min(corr)>scenario.corr_th
                     break;
                 elseif min(corr)>corr_tmp
                     corr_tmp = min(corr);
                     h_data_tmp = h_data;
                     h_pilot_tmp = h_pilot;
                     UE_index_tmp = UE_index;
                 else
                     h_data = h_data_tmp;
                     h_pilot = h_pilot_tmp;
                     UE_index = UE_index_tmp;
                 end
             end
             
             if scenario.corr_flag == 0
                 if max(corr)<scenario.corr_th
                     break;
                 elseif max(corr)<corr_tmp
                     corr_tmp = max(corr);
                     h_data_tmp = h_data;
                     h_pilot_tmp = h_pilot;
                     UE_index_tmp = UE_index;
                 else
                     h_data = h_data_tmp;
                     h_pilot = h_pilot_tmp;
                     UE_index = UE_index_tmp;
                 end
             end
             
        end
        
    end
else
    UE_index = scenario.precalculated_indices(scenario.index, scenario.seed,:);
    if ~all(UE_index)
        h_pilot_out = NaN; 
        h_data_out = NaN;
        h_srs_out = NaN;
        return
    end
    % directories were added at indices precalculation
    
    for j=1 : N_user
        
        if j<=N_target
            comb = scenario.comb;
        else
            comb = scenario.int_comb;
        end
        first_ind = 1;
        if comb && scenario.comb_split
            first_ind = 2-mod(j,2); % 1 for odd, 2 for even
        end
        index=num2str(UE_index(j));
        H_tmp = load_channel(index, Nrx, 0);        
        
        % Extract channel for pilot symbols
        h_pilot(j,:,first_ind:(1+comb):N_used,1) = (sqrt(1+comb))*H_tmp(1,1:Nrx,first_ind:(1+comb):N_used, pilot_positions(1));
        h_pilot(j,:,first_ind:(1+comb):N_used,2) = (sqrt(1+comb))*H_tmp(1,1:Nrx,first_ind:(1+comb):N_used, pilot_positions(2));
        
        % Extract channel for data symbols
        h_data(j,:,1:N_used, data_positions1-0) = H_tmp(1,1:Nrx,1:N_used, data_positions1);
        h_data(j,:,1:N_used, data_positions2-1) = H_tmp(1,1:Nrx,1:N_used, data_positions2);
        h_data(j,:,1:N_used, data_positions3-2) = H_tmp(1,1:Nrx,1:N_used, data_positions3);

    end 
end

comb_SRS=0;
N_used_SRS=16*RB_size;
h_srs=zeros(N_user,Nrx,N_used_SRS,N_pilot_sym);

for j=1:N_user
    
    index=num2str(UE_index(j));
    H_tmp = load_channel(index, Nrx, scenario.SRS_delay);
    % Extract srs symbols
    h_srs(j,:,1:(1+comb_SRS):N_used_SRS,1) = (sqrt(1+comb_SRS))*H_tmp(1,1:Nrx,1:(1+comb_SRS):N_used_SRS, pilot_positions(1));
    h_srs(j,:,1:(1+comb_SRS):N_used_SRS,2) = (sqrt(1+comb_SRS))*H_tmp(1,1:Nrx,1:(1+comb_SRS):N_used_SRS, pilot_positions(2));
    
end

% disp(UE_index);

h_srs_out=h_srs;
h_pilot_out=h_pilot;
h_data_out=h_data;
end

function H = load_channel(index, Nrx, SRS)
    %fast channel loading
    if 1%Nrx ~= 1024
        if SRS
            H = load(['SRS_chan' '_seed' index]).H_new;
        else
            H = load(['temp_chan' '_seed' index]).H_new;
        end
    else
        H = load(['data_50G_1024RX_LOS_seed' index]).H_new;
    end
    
    if numel(size(H)) < 4
        H = repmat(H, [1,1,1,14]);
    end
end