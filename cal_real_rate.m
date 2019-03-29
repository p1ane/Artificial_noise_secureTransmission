clc;
clear;
load Omega_freq_128
load Lambda

CCM = Omega_freq_128(:,:,:,1:8);
CCMev = Omega_freq_128(:,:,:,9);

[Nr,Nt,~,Nu] = size(CCM);
[Ne,~,~] = size(CCMev);

NumSamples = 100;
h_freq = reform_H_beam(CCM, Nt, Nr, NumSamples, Nu);
h_freq_eve = reform_Heve_beam(CCMev,Nt,Ne,NumSamples);

SNR = -10:5:20;
NSNR = length(SNR);
rate_Real = zeros(1,NSNR);
rate_de_withAN = zeros(1,NSNR);
for nSNR = 1:NSNR
    rate_de_withAN(nSNR) = cal_Ergodic_rate(Lambda_optimal_withAN(:,:,nSNR),CCM,CCMev,h_freq,h_freq_eve);
    rate_Real(nSNR) = cal_rate_real(Lambda_optimal_withAN(:,:,nSNR),CCM,CCMev,h_freq,h_freq_eve);
end