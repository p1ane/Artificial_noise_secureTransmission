clear all;
clc;

SNR = -10:5:20;    %–≈‘Î±»∑∂Œß
NSNR = length(SNR);

load result;
load rate_Real;

figure;
hold on;


plot(SNR,rate_Real/log(2),'b-','LineWidth',1.5);
plot(SNR,rate_Monte_withAN/log(2),'r--','LineWidth',1.5);
plot(SNR,rate_de_withAN/log(2),'ks','LineWidth',1.5);
plot(SNR,rate_de_withoutAN/log(2),'m-','LineWidth',1.5);

legend('With AN,Ergoric secrevy sum rate','With AN,Lower bound','With AN,DE of lower bound','Without AN','location','NorthWest');
box on;
grid minor;

hold off;

xlabel('SNR (dB)');
ylabel('Secrecy sum rate (bits/s/Hz)');

