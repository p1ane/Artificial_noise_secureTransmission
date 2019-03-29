clc;
clear;
load my_iter;

SNR = -10:10:20;
NSNR = length(SNR);
ITER = 0:9;
NITER = length(ITER);
figure;
hold on;



plot(ITER,rate_iter_withAN(4,:)/log(2),'m-^','LineWidth',1.5);
plot(ITER,rate_iter_withAN(3,:)/log(2),'k-v','LineWidth',1.5);
plot(ITER,rate_iter_withAN(2,:)/log(2),'b-o','LineWidth',1.5);
plot(ITER,rate_iter_withAN(1,:)/log(2),'r-s','LineWidth',1.5);

legend('20 dB','10 dB','0 dB','-10 dB','location','NorthWest');
box on;
grid minor;

xlabel('Number of iterations');
ylabel('Secrecy sum rate (bits/s/Hz)');

hold off;