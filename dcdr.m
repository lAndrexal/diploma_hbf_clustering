N_init = 154;
N_corrupted = 20;
in = randi([0 1],N_init,1);
coded = nrLDPCEncode(in,1);
p = randperm(numel(coded),N_corrupted); % positions of bits to corrupt
coded(p) = 1-coded(p);
coded_soft = double(1-2*coded);
[in_restored,numIters,checkBits] = nrLDPCDecode(coded_soft,1,25);
all(in_restored==in)