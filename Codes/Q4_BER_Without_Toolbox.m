% Simulating across a range of SNR
Eb_No = -3: 1: 10;
% Frame Length
bit_num = 10000;

% Main loop
for aa = 1: 1: length(Eb_No)
    % To convert Eb/No numbers to channel SNR, use the formula below.
    SNR = Eb_No + 10*log10(2);
    T_Errors = 0;
    T_bits = 0;
    
    
    % Continue until 100 errors.
    while T_Errors < 100
    
        % Make a few data bits
        uncoded_bits  = round(rand(1,bit_num));
        
        % For Quadrature Carriers, split the stream into two streams.
        B1 = uncoded_bits(1:2:end);
        B2 = uncoded_bits(2:2:end);
        
        % QPSK modulator set to pi/4 radians constellation
        % If you want to change the constellation angles
        % just change the angles.
        qpsk_sig = ((B1==0).*(B2==0)*(exp(i*pi/4))+(B1==0).*(B2==1)...
            *(exp(3*i*pi/4))+(B1==1).*(B2==1)*(exp(5*i*pi/4))...
            +(B1==1).*(B2==0)*(exp(7*i*pi/4))); 
        
        % Noise variance
        N0 = 1/10^(SNR(aa)/10);
        
        % To the recipient, send a Gaussian Link message.
        rx = qpsk_sig + sqrt(N0/2)*(randn(1,length(qpsk_sig))+i*randn(1,length(qpsk_sig)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % At the receiver, there is a QPSK demodulator.
        B4 = (real(rx)<0);
        B3 = (imag(rx)<0);
        
        uncoded_bits_rx = zeros(1,2*length(rx));
        uncoded_bits_rx(1:2:end) = B3;
        uncoded_bits_rx(2:2:end) = B4;
    
        % Bit Errors Calculation
        diff = uncoded_bits - uncoded_bits_rx;
        T_Errors = T_Errors + sum(abs(diff));
        T_bits = T_bits + length(uncoded_bits);
        
    end

 
    BER(aa) = T_Errors / T_bits;

end
  
% BER through Simulation
figure(1);
semilogy(SNR,BER,'or');
hold on;
xlabel('SNR in dB');
ylabel('BIT EROOR RATE');
title('In a Gaussian context, display SNR vs. BER for QPSK Modulation.');
% Theoretical BER
figure(1);
theoryBer = 0.5*erfc(sqrt(10.^(Eb_No/10)));
semilogy(SNR,theoryBer);
grid on;
legend('Simulated', 'Theoretical');


%%%%%%%%%%%%%%%% Bit Error Rate for 16-QAM modulation %%%%%%%%%%%%%%%%%%%

N = 10^5; % the total amount of symbols
M = 16;   % the size of the constellation
k = log2(M); 

% for 16-QAM
alphaRe = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
alphaIm = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
k_16QAM = 1/sqrt(10);

Eb_N0_dB  = [0:15]; % multiple Es/N0 values
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

% Conversion of gray codes
ref = [0:k-1];
map = bitxor(ref,floor(ref/2));
[tt ind] = sort(map);                                

for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    ipBit = rand(1,N*k,1)>0.5; 
    ipBitReshape = reshape(ipBit,k,N).';
    bin2DecMatrix = ones(N,1)*(2.^[(k/2-1):-1:0]) ; % Binary to decimal
    
    % real
    ipBitRe =  ipBitReshape(:,[1:k/2]);
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
    
    % imaginary
    ipBitIm =  ipBitReshape(:,[k/2+1:k]);
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    
    % constellations based on Gray coded symbols
    modRe = alphaRe(ipGrayDecRe+1);
    modIm = alphaIm(ipGrayDecIm+1);
    
    
    mod = modRe + j*modIm;
    s = k_16QAM*mod; 
    
    % noise
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; 
    
    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    
    % real
    y_re = real(y)/k_16QAM;
    % imaginary
    y_im = imag(y)/k_16QAM; 

    % rounding up to the next letter of the alphabet
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(find(ipHatRe>max(alphaRe))) = max(alphaRe);
    ipHatRe(find(ipHatRe<min(alphaRe))) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(find(ipHatIm>max(alphaIm))) = max(alphaIm);
    ipHatIm(find(ipHatIm<min(alphaIm))) = min(alphaIm);

    % Converting from Constellation to Decimal
    ipDecHatRe = ind(floor((ipHatRe+4)/2+1))-1; % LUT based
    ipDecHatIm = ind(floor((ipHatIm+4)/2+1))-1; % LUT based

    % converting a string to binary
    ipBinHatRe = dec2bin(ipDecHatRe,k/2);
    ipBinHatIm = dec2bin(ipDecHatIm,k/2);

    % translating a binary text to a numerical value
    ipBinHatRe = ipBinHatRe.';
    ipBinHatRe = ipBinHatRe(1:end).';
    ipBinHatRe = reshape(str2num(ipBinHatRe).',k/2,N).' ;
    
    ipBinHatIm = ipBinHatIm.';
    ipBinHatIm = ipBinHatIm(1:end).';
    ipBinHatIm = reshape(str2num(ipBinHatIm).',k/2,N).' ;

    % Counting both real and imaginary errors
    nBitErr(ii) = size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1) ;

end 
simBer = nBitErr/(N*k);
theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));

figure
semilogy(Eb_N0_dB,theoryBer,'bs-',Eb_N0_dB,simBer,'rx');
axis([0 15 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('BIT ERROR RATE')
title('For 16-QAM modulation, the bit error probability curve is shown.')