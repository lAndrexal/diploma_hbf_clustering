P = [ 36  39   39    0  -1  -1 
      44  37   -1   -1   0  -1    
       3  14   22   -1  -1   0  ];
blockSize = 48;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix)
cfgLDPCDec = ldpcDecoderConfig(pcmatrix)

M = 4;
maxnumiter = 10;
snr = [3 6 20];
numframes = 10;

ber = comm.ErrorRate;
ber2 = comm.ErrorRate;

for ii = 1:length(snr)
    for counter = 1:numframes
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1,'int8');
        % Transmit and receive with LDPC coding
        encodedData = ldpcEncode(data,cfgLDPCEnc);
        modSignal = pskmod(encodedData,M,InputType='bit');
        [rxsig, noisevar] = awgn(modSignal,snr(ii));
        demodSignal = pskdemod(rxsig,M, ...
            OutputType='approxllr', ...
            NoiseVariance=noisevar);
        rxbits = ldpcDecode(demodSignal,cfgLDPCDec,maxnumiter);
        errStats = ber(data,rxbits);
        % Transmit and receive with no LDPC coding
        noCoding = pskmod(data,M,InputType='bit');
        rxNoCoding = awgn(noCoding,snr(ii));
        rxBitsNoCoding = pskdemod(rxNoCoding,M,OutputType='bit');
        errStatsNoCoding = ber2(data,int8(rxBitsNoCoding));
    end
    fprintf(['SNR = %2d\n   Coded: Error rate = %1.2f, ' ...
        'Number of errors = %d\n'], ...
        snr(ii),errStats(1),errStats(2))
    fprintf(['Noncoded: Error rate = %1.2f, ' ...
        'Number of errors = %d\n'], ...
        errStatsNoCoding(1),errStatsNoCoding(2))
    reset(ber);
    reset(ber2);
end