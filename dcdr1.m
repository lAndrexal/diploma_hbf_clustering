% Parameters for LDPC
n = 288; % Length of codeword
k = 144; % Length of message
rate = k / n; % Code rate

% Step 1: Generate Parity-Check Matrix (H)
% H is (n-k) x n matrix
dv = 3; % Column weight
dc = (dv * n) / (n - k); % Row weight
H = make_ldpc(n, n - k, dc, dv);

% Step 2: Generate Generator Matrix (G)
[Ht, ~] = size(H);
P = H(:, 1:k);
I = eye(k); % Identity matrix
G = [I, P']; % Generator matrix

% Step 3: Encode Message
msg = randi([0, 1], 1, k); % Random message
codeword = mod(msg * G, 2); % LDPC encoding

% Step 4: Add Noise to Transmit Codeword
snr = 3; % Signal-to-noise ratio in dB
transmitted = 1 - 2 * codeword; % BPSK modulation
noisy = awgn(transmitted, snr, 'measured'); % Add Gaussian noise

% Step 5: Decode Received Codeword
L = 2 * noisy / (10^(snr / 10)); % LLR calculation for decoding
maxIter = 50; % Maximum number of iterations
decoded_msg = decode_ldpc(L, H, maxIter);

% Step 6: Evaluate Results
fprintf("Original message:     ");
disp(msg);
fprintf("Decoded message:      ");
disp(decoded_msg);
fprintf("Decoding successful:  %d\n", isequal(logical(msg), decoded_msg));

%% Supporting Functions

% Function to create LDPC parity-check matrix (H)
function H = make_ldpc(n, m, dc, dv)
    H = zeros(m, n);
    for i = 1:n
        row_indices = randperm(m, dv);
        H(row_indices, i) = 1;
    end
    % Ensure each row meets the weight requirement
    for j = 1:m
        while sum(H(j, :)) > dc
            col_indices = find(H(j, :) == 1);
            H(j, col_indices(end)) = 0;
        end
    end
end

% Function for LDPC decoding using Belief Propagation
function decoded = decode_ldpc(L, H, maxIter)
    [m, n] = size(H);
    R = zeros(m, n); % Initialize messages
    Lq = L; % Initialize intrinsic LLRs
    for iter = 1:maxIter
        % Update Check Nodes
        for i = 1:m
            idx = find(H(i, :) == 1);
            for j = idx
                R(i, j) = prod(tanh(Lq(idx(idx ~= j)) / 2));
            end
        end

        % Update Variable Nodes
        for j = 1:n
            idx = find(H(:, j) == 1);
            Lq(j) = L(j) + sum(R(idx, j));
        end

        % Hard Decision
        decoded = Lq < 0;

        % Check Parity
        if mod(H * decoded', 2) == 0
            break;
        end
    end
end
