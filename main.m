clc;
clear;


load Omega_freq_128
CCM = Omega_freq_128(:,:,:,1:8);
CCMev = Omega_freq_128(:,:,:,9);

[Nr,Nt,~,Nu] = size(CCM);
[Ne,~,~] = size(CCMev);

NumSamples = 200;
h_freq = reform_H_beam(CCM, Nt, Nr, NumSamples, Nu);
h_freq_eve = reform_Heve_beam(CCMev,Nt,Ne,NumSamples);

SNR = -10:5:20;
NSNR = length(SNR);
Lambda = zeros(Nt,Nu+1,NSNR);


rate_Monte_withAN = zeros(1,NSNR);
rate_de_withAN = zeros(1,NSNR);
rate_Monte_withoutAN = zeros(1,NSNR);
rate_de_withoutAN = zeros(1,NSNR);
Lambda_optimal_withAN = zeros(Nt,Nu+1,NSNR);
Lambda_optimal_withoutAN = zeros(Nt,Nu+1,NSNR);
NITER = 10;
rate_iter_withAN = zeros(NSNR,NITER);
rate_iter_withoutAN = zeros(NSNR,NITER);

for nSNR = 1:NSNR
    snr = 10^(SNR(nSNR)/10);
    Lambda(:,:,nSNR) = ones(Nt,Nu+1) / (Nu+1) /Nt  * snr;
    
 
    
    [Lambda_optimal_withAN(:,:,nSNR),rate_iter_withAN(nSNR,:)] = cal_CCCP_withAN(Lambda(:,:,nSNR),snr,CCM,CCMev,NITER);
    [Lambda_optimal_withoutAN(:,:,nSNR),rate_iter_withoutAN(nSNR,:)] = cal_CCCP_withoutAN(Lambda(:,:,nSNR),snr,CCM,CCMev,NITER);
   
    
    rate_de_withoutAN(nSNR) = cal_DE_rate(Lambda_optimal_withoutAN(:,:,nSNR),CCM,CCMev);
    rate_Monte_withoutAN(nSNR) = cal_Ergodic_rate(Lambda_optimal_withoutAN(:,:,nSNR),CCM,CCMev,h_freq,h_freq_eve);
    rate_de_withAN(nSNR) = cal_DE_rate(Lambda_optimal_withAN(:,:,nSNR),CCM,CCMev);
    rate_Monte_withAN(nSNR) = cal_Ergodic_rate(Lambda_optimal_withAN(:,:,nSNR),CCM,CCMev,h_freq,h_freq_eve);
end