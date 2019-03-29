clc;
clear;

load CCM1;
load CCMev1;

[Nr,Nt,Nu] = size(CCM);
[Ne,~] = size(CCMev);

NumSamples = 200;
h_freq = reform_H_beam(CCM, Nt, Nr, NumSamples, Nu);
h_freq_eve = reform_Heve_beam(CCMev,Nt,Ne,NumSamples);

SNR = 20;
NSNR = length(SNR);
Lambda = zeros(Nt,Nu+1,NSNR);


rate_Monte = zeros(1,NSNR);
rate_de = zeros(1,NSNR);


for nSNR = 1:NSNR
    snr = 10^(SNR(nSNR)/10);
    
    Lambda(:,:,nSNR) = ones(Nt,Nu+1) / (Nu+1) /Nt  * snr;      
    
end

for nSNR = 1:NSNR
    rate_de(nSNR) = cal_DE_rate(Lambda(:,:,nSNR),CCM,CCMev);
    rate_Monte(nSNR) = cal_Ergodic_rate(Lambda(:,:,nSNR),CCM,CCMev,h_freq,h_freq_eve);   
end