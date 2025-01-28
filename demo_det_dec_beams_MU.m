clear all; %#ok<CLALL> 

%%%%%%%%%%%%%%%%%%%%%%%%% Files parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix path for server
if isunix()
    decode_path = "~/Desktop/Research/SourceCodeOur/";
    channels_path = "~/Desktop/Research/Resources/";
else
    decode_path = "../../SourceCodeOur/";
    channels_path = "../../Resources/";
end
addpath(decode_path+"decoding/")
addpath(channels_path+"data_1024RX_rank1")
addpath("Beamforming")

%%%%%%%%%%%%%%%%%%%%%%%%% Define scenario %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters

scenario.SNR_range=-23:3:-8;    % SNR, dB
scenario.N_scenarios = 75; % channel realization 1...140
scenario.N_seeds = 2;      % seed 1...99 (ideally = 99 , but takes too much time)
ML_coef=0;                 

% The first UE has always a power of 0 dB !!!
scenario.UE_power=zeros(1,4);%[0];                 % power ratio between users, dB
scenario.N_target_UE=4;                         % number of target users, dB

% BS geometrical parameters

scenario.Nrx=1024;            % 1024 / 64 / 8 / 4 antennas

% Signal parameters

scenario.RB_num=4;         % number of resource blocks (frequency band)
scenario.RB_size=12;
scenario.N_TTI=1;          % number of TTIs (1 TTI = 14 symbols = 12 data + 2 pilot)
scenario.N_data_sym = 12;

scenario.QAM_order=4;
scenario.code_block_index=0; 
scenario.SC_FDMA=0; 

scenario.ZC_root=1:20;%[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];          % Zadoff-Chu root (high impact)
scenario.ZC_shift=zeros(1,24);          % Zadoff-Chu shift (low impact)

scenario.comb=1;            % comb location (0 - off, 1 - on)
scenario.int_comb=0;        % comb location for interference (0 - off, 1 - on)

% Parameters for user split

scenario.ideal_user_split = 1;

scenario.comb_split=1; % Allocate different users on different subcarriers

N_users = length(scenario.UE_power);
scenario.OC_code = ones(N_users,2);
tmp = 1:N_users; tmp = (mod(tmp,4) == 3) + (mod(tmp,4) == 0);
scenario.OC_code(find(tmp),2) = -1; % Orthogonal Cover Code 
scenario.OC_code(scenario.N_target_UE+1:end,:) = 1; % no OCC for interf

% Channels parameters

scenario.corr_th=0.8;     % correlation threshold
scenario.N_iter=1;        % max N-of-iter to achieve the desired correlation
scenario.corr_flag=1;     % 1  -   high correlation,   0  -  low correlation, 0<x<1 - no care

%%%%%%%%%%%%%%%%%%%%%%%%% Algorithms parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%

% CE parameters

scenario.ideal_CE=1; 
scenario.consider_interf_on_CE = 0;   
PHY_param.max_N_peaks=8;       % max N-of-peaks to compensate (4...12)   
PHY_param.N_stop = 8;         % max N-of-peaks to compensate (only for _PIC version)  
scenario.interp_order=0;

estimator_osin = @(IN_DATA, estimationParams, ML_coef) ...
    CE_TTI_det2_synch_beams(IN_DATA, estimationParams, ML_coef, PHY_param); % Skoltech CE algorithm
estimator_osin_spread = @(IN_DATA, estimationParams, ML_coef) ...
    CE_TTI_det2_beams(IN_DATA, estimationParams, ML_coef, PHY_param); % asynchronous version of Skoltech algorithm (only correlation with the first beam is used)
estimator_SRS = @(h_f_noisy, SRS_Params, ML_coef) ...
    CE_TTI_det2_SRS(h_f_noisy, SRS_Params, ML_coef, PHY_param); % Define beams basis estimation algorithm

% SRS beamspace parameters

scenario.beam_transform=1;
scenario.N_ports=32;
scenario.ideal_SRS = 1;
scenario.SRS_delay = 0; % Use SRS or DMRS for beamspace estimation
scenario.beamspace_mode = 2; % 2 - Kronecker product of DFTs
scenario.hybrid = 1;    % For 1024 antennas
scenario.analog_mode=1;
scenario.Ndigital = 64; % N of ADCs
scenario.new_beamspace = 0;
scenario.training_mode = 3;
scenario.digital_pol = 1;
scenario.subarray_phase = 0;

scenario.training = false;
% General detector parameters

scenario.RB_num_Ruu=4;         % number of resource blocks to estimate Ruu
scenario.ideal_Ruu=0;
scenario.interp_step=1;

% Half precision parameters (double parameters are in ML_coef var)

scenario.half_precision = 0;

scenario.alpha = 1;
scenario.singleRUU = 0;
scenario.col2Diag = 0; % 
scenario.stop_sorting = 1;
scenario.modeA = 0; % Cholesky (0) or LQ (-1) for H pseudoinverse
scenario.Wmode = 0; % Bitmask for partial single (0) - fully half
                    %                            (32) - fully single
scenario.col2LQ = 4;
scenario.single_sRuu=1;
scenario.single_sChol=1;
scenario.single_sQ=0;

%%%%%%%%%%%%%%%%%%%%%%%%% Define channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate channel indices for the simulation
if length(scenario.UE_power)>1 && scenario.N_iter>1
    scenario.precalculated_indices = precalculate_indices(scenario);
end
scenario.precalculated_indices = load("HC_sub_pi.mat").precalculated_indices;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% VALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_scenario = scenario; % To backup scenario after adjustment
description = {}; % What each calculated curve means in legend and save_file
curve_id = 0;

curve_id = curve_id+1;
description{curve_id} = "-4 dft 32"; % Use this as comment to FER call
scenario.analog_mode=-4;
scenario.beamspace_mode=2;
scenario.N_ports=32;
FER_skoltech{curve_id} = scenarios_validation(scenario, estimator_osin, estimator_SRS, ML_coef);
scenario = default_scenario;

curve_id = curve_id+1;
scenario.precalculated_indices = load("HC_sub_pi.mat").precalculated_indices;
scenario.clustering_path = "HC_sub";
description{curve_id} = "Clustering DFT"; % Use this as comment to FER call
scenario.analog_mode=-7;
scenario.beamspace_mode=-1;
scenario.N_ports=4;
FER_skoltech{curve_id} = scenarios_validation(scenario, estimator_osin, estimator_SRS, ML_coef);
scenario = default_scenario;
 
% curve_id = curve_id+1;
% description{curve_id} = "antenna";
% scenario.beam_transform=0;
% FER_skoltech{curve_id} = scenarios_validation(scenario, estimator_osin, estimator_SRS, ML_coef);
% scenario = default_scenario;


save("paper_reconstruct")
save_mat = 0;
if save_mat; save(save_filename); end

%% %%%%%%%%%%%%%%%%%%% PLOT target UEs performance %%%%%%%%%%%%%%%%%%%%%%%%

% Define styles for curves
descriptors = ["r-*",... % 1
               "k-",...  % 2              
               "b-",...  % 3            
               "k-*",... % 4
               "m-",...  % 5
               "r-",...  % 6
               ];
%
SNR_range = scenario.SNR_range;
if length(SNR_range) > 1
    for i=1:scenario.N_target_UE
        figure(i);
        % figure(i,'Position',[0 0 500 500]);
        
        p = [];
        for j = 1:numel(FER_skoltech)
            
            if scenario.N_target_UE>1
                FER = FER_skoltech{j}(i,:);
            else
                FER = FER_skoltech{j};
            end
            
            p_ = semilogy(SNR_range, FER, descriptors(j),'LineWidth',2);
            p = [p p_];
            if j == 1
                hold on; grid on; grid minor; set(gca,'fontsize', 14); xlabel('SNR, dB'); xticks(SNR_range); ylabel('FER');
            end
        end
        
        axis([SNR_range(1) SNR_range(length(SNR_range)) 10^-2 1]);

        switch scenario.code_block_index; case 2 ; code_name=', LDPC(216,288)'; case 1 ; code_name=', LDPC(1152,2304)'; case 0 ; code_name=', LDPC(144,288)'; end

        UEs_info = num2str(scenario.N_target_UE) + " targets, "+num2str(numel(scenario.UE_power) - scenario.N_target_UE)+" interference, I/N=30dB, ";

        title(UEs_info+"User"+num2str(i)+", "+num2str(scenario.Nrx)+"RX, "+num2str(scenario.N_ports)+"beams (real SRS estimator), "+num2str(scenario.RB_num)+"RB, QAM"+num2str(2^scenario.QAM_order) + code_name);
        
        legend(p, description)
        
        legend('Location', 'southwest');
    end
end



% function for LEARNING 
function sum_error = scenarios_learning(scenario, estimator_DMRS, estimator_SRS, ML_coef, range) %#ok<DEFNU> 
k=0;

N_scenarios=scenario.N_scenarios;
N_seeds=scenario.N_seeds;
N_target_UE=scenario.N_target_UE;

scen=scenario;
SNR_error=0;

for SNR=range
    
    scen.SNR=SNR;
    err_seed=zeros(N_target_UE,N_scenarios*N_seeds);
    

    for j=0:(N_scenarios*N_seeds)-1
        scen.seed=rem(j,N_seeds)+1;
        scen.index=fix(j/N_seeds)+1;
        s(j+1)=scen; %#ok<AGROW> 
    end
    
    for j=1:(N_scenarios*N_seeds)
        sys=struct(s(j));
        [~,err_seed(:,j)]=tester_det_dec_beams_SCFDMA_MU(sys, estimator_DMRS, estimator_SRS, ML_coef);
    end
    
    err_scen=sum(err_seed,2)/N_seeds;
    
    SNR_error=SNR_error+(1.5^k)*err_scen;
    
    k=k+1;
    
end

sum_error=SNR_error;

end

% function for VALIDATION 
function sum_error = scenarios_validation(scenario, estimator_DMRS, estimator_SRS, ML_coef)
    
    N_scenarios=scenario.N_scenarios;
    N_seeds=scenario.N_seeds;
    N_target_UE=scenario.N_target_UE;
    
    scen=scenario;
    
    n_SNR = numel(scen.SNR_range);
    err_seed=zeros(N_scenarios*N_seeds, N_target_UE, n_SNR);
    err1=zeros(N_scenarios*N_seeds, N_target_UE, n_SNR);
    err2=zeros(N_scenarios*N_seeds, N_target_UE, n_SNR);
    for j=0:(N_scenarios*N_seeds)-1
        scen.seed=rem(j,N_seeds)+1;
        scen.index = fix(j/N_seeds)+1;
        s(j+1)=scen; %#ok<AGROW> 
    end
    ind_scen = 0;
    parfor j=1:(N_scenarios*N_seeds)%
        sys=struct(s(j));
        sys.N_user=length(sys.UE_power);       % number of all users in scenario
        
        seed = 100*sys.seed+10000*sys.index;
        rng(seed);

        channels = struct;
        [channels.h_pilot, channels.h_data, channels.h_srs] = get_channel(sys); % read channel
        if isnan(channels.h_pilot)
            for i = 1:n_SNR % Do not vectorize, otherwise parfor error
                err1(j,:,i) = NaN;
                err2(j,:,i) = NaN;
                err_seed(j,:,i) = NaN;
            end
            continue
        end
        ind_scen = ind_scen + 1;
        
        SNR_range = sys.SNR_range;
        for i = 1:n_SNR
            sys.SNR = SNR_range(i);
            [err1(j,:,i),err2(j,:,i),err_seed(j,:,i)]=tester_det_dec_beams_SCFDMA_MU(sys, estimator_DMRS, estimator_SRS, ML_coef, channels);
        end
    end
    
    err_scen=squeeze(mean(err_seed,1,'omitnan'));
    err1 = squeeze(mean(err1,1,'omitnan'));
    err2 = squeeze(mean(err2,1,'omitnan'));

    if scenario.half_precision
        disp(["Mode: ", scenario.half_mode, ", ideal CE: ", scenario.ideal_CE]);
    else
        disp(["Mode: single", ", ideal CE: ", scenario.ideal_CE]);
    end
    disp("MSE channel:")
    disp(err1);
    disp("MSE detector:")
    disp(err2);
    disp("FER:")
    disp(err_scen);
    
    sum_error = err_scen;
end



