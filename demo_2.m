SNR = 0:2:20;
BER = zeros(11,1);
SER = zeros(11,1);
for i=1:11
    [SER(i,1) , BER(i,1) , out] = mpam(2^18,2^3,SNR(i)); 
end

scatter((1:1:length(BER)), BER);
