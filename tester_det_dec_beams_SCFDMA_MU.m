function [RMSE_CE, MSE_detector, FER] = tester_det_dec_beams_SCFDMA_MU(scenario, estimator_DMRS, estimator_SRS, ML_coef1, channels)

% Same as in demo to work in training mode
N_user=length(scenario.UE_power);       % number of all users in scenario
scenario.N_user = N_user;           

% If channels are not passed to function - load
if ~exist("channels",'var') || ~isfield(channels,'h_srs')
    seed = 100*scenario.seed+10000*scenario.index;
    rng(seed);
    [h_pilot, h_data, h_srs,~]  = get_channel(scenario);    % read channel
else 
    h_pilot = channels.h_pilot; 
    h_data = channels.h_data;
    h_srs = channels.h_srs;
end

Nrx = scenario.Nrx; % N-of-Rx antennas
comb = scenario.comb;
SC_FDMA = scenario.SC_FDMA;

% RB_size = 12 subcarriers
RB_size=scenario.RB_size;
RB_num=scenario.RB_num;
RB_num_Ruu=scenario.RB_num_Ruu;

% N-of-subcarriers (max=600)
N_used=RB_num*RB_size; % var - use it to change simulation set 

% N of symbols per TTI
% N_ofdm = 14;

N_data_sym=scenario.N_data_sym;
% N-of-pilot-symbols
N_pilot_sym = 2;

% Pilot positions in TTI
pilot_positions=[4 11];

% Constelation order
QAM_order=scenario.QAM_order;
QAM_points=2^QAM_order;
SNR = scenario.SNR;


N_target_UE  = scenario.N_target_UE;    % number of target users in scenario
% N_scenarios=scenario.N_scenarios;

ZC_root=scenario.ZC_root;
ZC_shift=scenario.ZC_shift;

% initialize LDPC
% switch scenario.code_block_index
%     case 2
%         matrix_name='LDPC_matrix_72_288';
%         code_rate=0.75;
%         code_length=288;
%     case 1
%         matrix_name='LDPC_matrix_1152_2304';
%         code_rate=0.5;
%         code_length=2304;
%     case 0
%         matrix_name='LDPC_matrix_144_288';
%         code_rate=0.5;
%         code_length=288;
% end

K_LDPC = 144;
code_length=288;
code_rate = K_LDPC/code_length;
P = [ 36  39   39    0  -1  -1 
      44  37   -1   -1   0  -1    
       3  14   22   -1  -1   0  ];
blockSize = 48;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);

% HS = load(matrix_name).HS;
% hLDPCenc = comm.LDPCEncoder(HS);
% init_ldpc = @(x) decode_soft(0, x);
% [H_LDPC, ~] = alist2sparse([matrix_name '.alist']);
% [temp, N_LDPC] = size(H_LDPC);
% K_LDPC = N_LDPC - temp;
% [ldpc, ~, ~] = init_ldpc([matrix_name '.alist']);
% scale_array = 0.75*ones(1, N_LDPC);
% offset_array = 0*ones(1,N_LDPC);


% calculate default averaged power (from the file)
UE_power_default = zeros(1,N_user);
for i=1:N_user
    UE_power_default(i)=mean( mean( mean( squeeze(h_pilot(i,:,:,:)).*conj(squeeze(h_pilot(i,:,:,:))) )));
end

% generate white noise for pilots
noise_p = (randn(Nrx, N_used, N_pilot_sym)+1i*randn(Nrx, N_used, N_pilot_sym)) / sqrt(2);
% generate white noise for data symbols
noise_d = (randn(Nrx, N_used, N_data_sym)+1i*randn(Nrx, N_used, N_data_sym)) / sqrt(2);
% estimate noise power per subcarrier
noise_power = sum(sum(sum(abs(noise_p).*abs(noise_p))))/(Nrx*N_used*N_pilot_sym); %   -  algorithm is required !!!!!


h_pilot_noisy=noise_p;
h_pilot_noisy_splitted(1:N_user, :, :, :) = repmat(permute(noise_p,[4,1,2,3]), [N_user, 1,1,1]); 
h_data_noisy=noise_d;
v=0 : N_used-1;
data_uncoded = zeros(N_user, N_used*QAM_order*N_data_sym / code_length, K_LDPC);
data_coded = zeros(N_used*QAM_order*N_data_sym / code_length, code_length);
s_tx_f = zeros(N_user, N_used, N_data_sym);

for i=1:N_user
    
    ZC_vec = circshift( exp(-1i*(pi*ZC_root(i)*v.*(v+1)) / N_used)  ,  [0 ZC_shift(i)]   );
    ZC_array = repmat(ZC_vec,Nrx,N_pilot_sym); 
    ZC_array = reshape(ZC_array,Nrx,N_used,2);
    
    if i<=N_target_UE
        tmp_SNR=SNR;
    else
        tmp_SNR=0; % ATTENTION! CHANGED  FROM "0"!!!!
    end
    gain = sqrt(10^(  (scenario.UE_power(i)+tmp_SNR)  /  10))  / sqrt(UE_power_default(i));
    
    h_pilot(i,:,:,:) = gain*squeeze(h_pilot(i,:,:,:));
    h_data(i,:,:,:) = gain*squeeze(h_data(i,:,:,:));
    
    for j=1:N_used*QAM_order*N_data_sym / code_length
        data_uncoded(i,j,:) = randi([0 1], K_LDPC, 1);
        tmp=squeeze(data_uncoded(i,j,:));
        data_coded(j,:) = ldpcEncode(tmp,cfgLDPCEnc);
    end
    
    % Modulate data
    dataInMatrix = reshape(data_coded,[],QAM_order);                 % Reshape data into binary k-tuples, k = log2(M)
    dataSymbolsIn = bi2de(dataInMatrix);                             % Convert to integers
    tmp=qammod(dataSymbolsIn.',QAM_points,'UnitAveragePower', true); % Gray coding, phase offset = 0;
    s_tx = reshape( tmp , N_used, [] );

    % generate channel response of single UE with noise
    data_mod = zeros(Nrx, N_used, N_data_sym);
    for k=1:N_data_sym
        
        if SC_FDMA
            s_tx_f(i,:,k)=fft(s_tx(:,k))/sqrt(N_used);
        else
            s_tx_f(i,:,k)=s_tx(:,k);
        end
        
        
        for m=1:N_used
            data_mod(:,m,k)=squeeze(h_data(i,:,m,k))*s_tx_f(i,m,k);
        end
    end
    
    h_pilot_noisy=h_pilot_noisy+squeeze(h_pilot(i,:,:,:)).*ZC_array.* reshape(scenario.OC_code(i,:), [1,1,N_pilot_sym]);
    h_pilot_noisy_splitted(i,:,:,:) = squeeze(h_pilot_noisy_splitted(i,:,:,:)) + ...
        squeeze(h_pilot(i,:,:,:)).*ZC_array.* reshape(scenario.OC_code(i,:), [1,1,N_pilot_sym]);
    % test interf impact on CE 
    if i == N_target_UE
        h_pilot_noisy_target_UE_only = h_pilot_noisy;
    end

    h_data_noisy=h_data_noisy+data_mod;
    
end

% dummy scaling for CE unit (depends on pilots scaling in comb mode)
h_pilot_noisy=h_pilot_noisy/sqrt(1+comb);
h_pilot_noisy_target_UE_only=h_pilot_noisy_target_UE_only/sqrt(1+comb);
h_pilot_noisy_splitted = h_pilot_noisy_splitted / sqrt(1+comb);
%LS CE calculation for DMRS
h_pilot_LS = zeros(N_target_UE, Nrx, N_used, N_pilot_sym);
for i=1:N_target_UE
    
    ZC_vec = circshift( exp(-1i*(pi*ZC_root(i)*v.*(v+1)) / N_used)  ,  [0 ZC_shift(i)] );
    ZC_array = repmat(ZC_vec,Nrx,N_pilot_sym); 
    ZC_array = reshape(ZC_array,Nrx,N_used,2);
    
    if scenario.consider_interf_on_CE == 1
        noisy_pilot = h_pilot_noisy;
    else
        noisy_pilot = h_pilot_noisy_target_UE_only;
    end
    if scenario.ideal_user_split
        h_pilot_LS(i,:,:,:)=squeeze(h_pilot_noisy_splitted(i,:,:,:)).*conj(ZC_array).* reshape(scenario.OC_code(i,:), [1,1,N_pilot_sym]);
    else
        h_pilot_LS(i,:,:,:)=noisy_pilot.*conj(ZC_array).* reshape(scenario.OC_code(i,:), [1,1,N_pilot_sym]);
    end
end

% DMRS params
DMRS_Params.SNR_dummy=SNR;
DMRS_Params.RB_size = RB_size;
DMRS_Params.RB_num = RB_num;
DMRS_Params.N_pilot = N_pilot_sym;
DMRS_Params.pilot_positions=pilot_positions;
DMRS_Params.comb=comb;
DMRS_Params.Nrx = scenario.Nrx;
DMRS_Params.N_ports = scenario.N_ports;
DMRS_Params.beam_transform = scenario.beam_transform;
DMRS_Params.comb_split = scenario.comb_split;
DMRS_Params.OC_code = scenario.OC_code;

DMRS_Params.interp_order = scenario.interp_order;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SRS params
SRS_Params=DMRS_Params;
SRS_Params.N_ports = scenario.N_ports;
SRS_Params.RB_num=16;  
N_ports=scenario.N_ports;

gain_SRS=sqrt(2);  % SRS power is higher then DMRS
if ~scenario.ideal_SRS

    comb_SRS=0;
    
    N_used_SRS=SRS_Params.RB_num*RB_size;

    % generate white noise for SRS
    h_srs_noisy = (randn(Nrx, N_used_SRS, N_pilot_sym)+1i*randn(Nrx, N_used_SRS, N_pilot_sym)) / sqrt(2);
    v=0 : N_used_SRS-1;
    for i=1:N_user
        
        ZC_vec = circshift( exp(-1i*(pi*ZC_root(i)*v.*(v+1)) / N_used_SRS)  ,  [0 ZC_shift(i)]   );
        ZC_array = repmat(ZC_vec,Nrx,N_pilot_sym); 
        ZC_array = reshape(ZC_array,Nrx,N_used_SRS,2);
        
        if i<=N_target_UE
            tmp_SNR=SNR;
        else
            tmp_SNR=0;
        end
        gain = sqrt(10^(  (scenario.UE_power(i)+tmp_SNR)  /  10))  / sqrt(UE_power_default(i));
        
        h_srs(i,:,:,:) = gain_SRS*gain*squeeze(h_srs(i,:,:,:));
        h_srs_noisy=h_srs_noisy+squeeze(h_srs(i,:,:,:)).*ZC_array;
        
    end
    
    % dummy (depends on pilots scaling in comb mode)
    h_srs_noisy=h_srs_noisy/sqrt(1+comb_SRS);

    % LS CE for calculation SRS
    h_srs_MU = zeros(N_user, Nrx, N_used_SRS);
    for i=1:N_user
        ZC_vec = circshift( exp(-1i*(pi*ZC_root(i)*v.*(v+1)) / N_used_SRS)  ,  [0 ZC_shift(i)]  );
        ZC_array = repmat(ZC_vec,Nrx,N_pilot_sym); 
        ZC_array = reshape(ZC_array,Nrx,N_used_SRS,2);
        h_srs_LS = h_srs_noisy.* conj(ZC_array);
        
        % Beam angles estimation via SRS (100 TTIs delay between SRS and current DMRS)
        
        [~,~,h_srs_f]=estimator_SRS(h_srs_LS,SRS_Params,ML_coef1);
        h_srs_MU(i,:,:)=h_srs_f;
        
    end
    [SRS_transform_matrix]=beamspace_selection(h_srs_MU, scenario);
else
    
    gain = sqrt(10.^(  (scenario.UE_power+tmp_SNR)  /  10))  ./ sqrt(UE_power_default);
    h_srs = gain_SRS.*reshape(gain, [N_user,1,1,1]).*h_srs;
    h_srs = mean(h_srs, 4);
    
    % default beamspace selection algorithm has power equalization
    if 0
        [SRS_transform_matrix]=beamspace_selection_new(h_srs, scenario);
    else
        [SRS_transform_matrix]=beamspace_selection(h_srs, scenario);
    end
end

% Channel transfer to the beam domain
h_beam_noisy = zeros(N_used, N_ports, N_pilot_sym);
h_beam_LS = zeros(N_target_UE, N_used, N_ports, N_pilot_sym);
for j=1:N_pilot_sym
    h_beam_noisy(:,:,j)=squeeze(h_pilot_noisy(:,:,j)).'*conj(SRS_transform_matrix);
    for i=1:N_target_UE
        % 64 antennas -> N_ports
        h_beam_LS(i,:,:,j)=squeeze(h_pilot_LS(i,:,:,j)).'*conj(SRS_transform_matrix);
    end
end

% specify input data for Channel Estimation unit

IN_DATA.SRS_transform_matrix=SRS_transform_matrix;
scenario.SRS_transform_matrix=SRS_transform_matrix;
transform_matrix=SRS_transform_matrix;

% Initialize both to pass into some functions
beam_amplitudes = zeros(N_target_UE, N_used, N_ports, N_data_sym);
h_data_recovered_f = zeros(N_target_UE, Nrx, N_used, N_data_sym);

RMSE_CE(1:N_target_UE) = 0;
if scenario.ideal_CE
    % ideal CE
    h_data_recovered_f = h_data(1:N_target_UE,:,:,:);
    for j=1:N_data_sym
        for i=1:N_target_UE
            % 64 antennas -> N_ports
            beam_amplitudes(i,:,:,j)=squeeze(h_data(i,:,:,j)).'*conj(SRS_transform_matrix);
        end
    end
    %transform_matrix=eye(Nrx);
    SNR_CE=Nrx*RB_num*(10^(SNR/10));
    
else
    beam_amplitudes_ideal = zeros(size(beam_amplitudes));
    for i=1:N_target_UE
        % channel estimation

        IN_DATA.h_f_noisy=squeeze(h_pilot_LS(i,:,:,:));
        IN_DATA.h_beam_noisy=squeeze(h_beam_LS(i,:,:,:));
        IN_DATA.user_ind = i;
        
        CE_DATA = estimator_DMRS(IN_DATA, DMRS_Params, ML_coef1);
        
        if i==1
            SNR_CE=CE_DATA.SNR; 
        end
        
        switch scenario.beam_transform
            case 1
                beam_amplitudes(i,:,:,:)=CE_DATA.SRS_beam_amplitudes;
                for j = 1:N_data_sym
                    beam_amplitudes_ideal(1,:,:,j)=squeeze(h_data(i,:,:,j)).'*conj(SRS_transform_matrix);
                end
                RMSE_CE(i) = sum(abs(beam_amplitudes(i,:,:,:) - beam_amplitudes_ideal(1,:,:,:)).^2,'all') / sum(abs(beam_amplitudes_ideal).^2,'all');
            case 0
                h_data_recovered_f(i,:,:,:)=CE_DATA.h_data_recovered_f;
                RMSE_CE(i) = sum(abs(h_data_recovered_f(i,:,:,:) - h_data(1:N_target_UE,:,:,:)).^2,'all') / sum(abs(h_data(1:N_target_UE,:,:,:)).^2,'all');
        end
    end
    
end

% MMSE Weight vector calculation
err_data=zeros(1,N_target_UE);
% error calculation

switch scenario.beam_transform
    case 1
        dim=N_ports;
        pilot_noisy=permute(h_beam_noisy, [2, 1, 3]);
        pilot_recovered=permute(beam_amplitudes, [1, 3, 2, 4]);
    case 0
        dim=Nrx;
        pilot_noisy=h_pilot_noisy;
        pilot_recovered=h_data_recovered_f;
end

pilot_recovered = pilot_recovered(:,:,:,1:N_pilot_sym);

pilot_noisy = pilot_noisy * sqrt (1+comb);

if comb
    for i=1:N_target_UE
        comb_matrix=ones(1, dim,N_used,N_pilot_sym);
        if scenario.comb_split
            first_ind = mod(i,2)+1;
        else
            first_ind = 2;
        end
        comb_matrix(1,:,first_ind:2:N_used,:)=zeros(1,dim,N_used/2,N_pilot_sym);
        pilot_recovered(i,:,:,:)=sqrt(2)*pilot_recovered(i,:,:,:).*comb_matrix;
    end
end

v=0 : N_used-1;
% s_rx_f = zeros(N_target_UE, RB_num_Ruu*RB_size, N_data_sym);
s_rx_f = zeros(N_target_UE, N_used, N_data_sym);
m=0;
for Ruu_index=1:round( N_used  /  (RB_num_Ruu*RB_size)  )
   
   U = zeros(RB_num_Ruu*RB_size, N_pilot_sym, dim);
   
   Ruu=zeros(dim,dim); 
   
    for i=1:N_pilot_sym
        for j=1:RB_num_Ruu*RB_size
            m=m+1;
            sc_index = RB_num_Ruu*RB_size*(Ruu_index-1)+j;
            
            if scenario.ideal_Ruu == 1
                
                if N_target_UE == N_user  % if there is no interference
                    Ruu=noise_power*eye(dim);
                else
                    
                    Ruu=Ruu+noise_power*eye(dim);
                    u= zeros(dim,1);
                    for k=N_target_UE+1:N_user
                        ZC = circshift( exp(-1i*(pi*ZC_root(k)*v.*(v+1)) / N_used)  ,  [0 ZC_shift(k)]   );
                        CE = squeeze(h_pilot(k,:,sc_index,i)).';
                        if scenario.beam_transform == 1
                            CE = SRS_transform_matrix' * CE;
                        end
                        u=u + CE * ZC(sc_index)*scenario.OC_code(k,i);
                    end
                    
                end
                
            else
                
                u= squeeze(pilot_noisy(:,sc_index,i));
                for k=1:N_target_UE
                    ZC = circshift( exp(-1i*(pi*ZC_root(k)*v.*(v+1)) / N_used)  ,  [0 ZC_shift(k)]   );
                    CE = squeeze(pilot_recovered(k,:,sc_index,i)).';
                    u=u - CE * ZC(sc_index)*scenario.OC_code(k,i);
                end
                
            end
            U(j,i,:) = u;
            Ruu=Ruu+u*u';
                            
        end
    end
    U = reshape(U, [RB_num_Ruu*RB_size*N_pilot_sym, dim]).';
    
    alpha = single(scenario.alpha) * single(noise_power); %ATTENTION
    U = single(U);
    if scenario.alpha ~= 0
        Ruu = U * U' / ( N_pilot_sym*RB_num_Ruu*RB_size ) + alpha * eye(dim);
    else
        Ruu = diag(repmat(mean(diag(U * U' / ( N_pilot_sym*RB_num_Ruu*RB_size ) )),[dim 1]));
    end
    
    for k=1:N_data_sym
        w = zeros(RB_num_Ruu*RB_size, N_target_UE, dim);
        for j=1:RB_num_Ruu*RB_size
            m=RB_num_Ruu*RB_size*(Ruu_index-1)+j;
            [~,H] = get_data_for_detector(m, k, beam_amplitudes, h_data_recovered_f, h_data_noisy, transform_matrix, scenario);
            x = [1:scenario.interp_step:(RB_num_Ruu*RB_size-1) RB_num_Ruu*RB_size];
            if ismember(j,x)
                R = chol(Ruu);
                H = R'\H;
                
                R2 = chol(H'*H + single(eye(N_target_UE)));
                w(j, :, :) = R2\(R2'\(H'))/R';
            end
        end
        if scenario.interp_step > 1
            for i = 1 : N_target_UE
                for j = 1 : dim
                    w(:, i, j) = interp1(x, w(x, i, j), 1:RB_num_Ruu*RB_size);
                end
            end
        end
        % Code for single precision or double precision detector
        for j=1:RB_num_Ruu*RB_size
            m=RB_num_Ruu*RB_size*(Ruu_index-1)+j;
            [Y,~] = get_data_for_detector(m, k, beam_amplitudes, h_data_recovered_f, h_data_noisy, transform_matrix, scenario);
            
            w_eq = reshape(w(j, :, :), [N_target_UE, dim]);

            s_rx_f(:,m,k) = w_eq*Y;
        end
    end
        
end

err_data=sum(abs(s_rx_f-s_tx_f(1:N_target_UE, :,:)).^2, [2,3]);
MSE_detector=err_data/(N_used*N_data_sym);
SNR_CE=max(10^-5, SNR_CE);


BER = ones(1,N_target_UE);
FER = ones(1,N_target_UE);
for i=1:N_target_UE

    if SC_FDMA
        s_rx=ifft(squeeze(s_rx_f(i,:,:)))*sqrt(N_used);
    else
        s_rx=squeeze(s_rx_f(i,:,:));
    end
    
    LLR_points = qamdemod(s_rx,QAM_points,'UnitAveragePower',true,'OutputType', 'approxllr', 'NoiseVariance', 1/SNR_CE);
    LLR_array = flip(permute( reshape (LLR_points, QAM_order, 1,[]), [3,2,1]), 3);
    
    LLR_vec = reshape( LLR_array, [], code_length);
    
    frame_err=0;
    data_decoded = zeros(QAM_order*N_data_sym*N_used/code_length, K_LDPC);
    for j=1:QAM_order*N_data_sym*N_used/code_length

        data_decoded(j,:) =  ldpcDecode(squeeze(LLR_vec(j,:).'),cfgLDPCDec,20);
        frame_err=frame_err+(sum(squeeze(data_decoded(j,:))~=squeeze(data_uncoded(i,j,:)).')>0);
    end
    BER(i) = sum(sum(data_decoded~=squeeze(data_uncoded(i,:,:))))/(N_used*N_data_sym*QAM_order*code_rate);
    FER(i)= frame_err/(QAM_order*N_data_sym*N_used/code_length);

end