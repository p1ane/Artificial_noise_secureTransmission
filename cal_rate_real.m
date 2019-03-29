function [sum_rate] = cal_rate_real(lambda,Omega,Omega_eve,h_freq,h_freq_eve)
[Nr,Nt,Nu] = size(Omega);
[Ne,~] = size(Omega_eve);
[K] = cal_K(lambda,Omega);

[~,~,~,NumSamples,~] = size(h_freq);
[~,~,~,NumSamples_eve] = size(h_freq_eve);
rate_all = zeros(Nu,1);

for k = 1:Nu
    
    rate_k_para1 = 0;
    for n_sample = 1:NumSamples
        rate_k_para1 = rate_k_para1 + log(abs(det( K(:,:,k) + h_freq(:,:,1,n_sample,k) * diag( lambda(:,k)+lambda(:,Nu+1)) * h_freq(:,:,1,n_sample,k)')));
        
    end
    rate_k_para1 = rate_k_para1 / NumSamples;
    
    rate_k_para2 = 0;
    for n_sample = 1:NumSamples_eve
        rate_k_para2 = rate_k_para2 + log(abs(det( eye(Ne) + h_freq_eve(:,:,1,n_sample) * diag(lambda(:,Nu+1)) * h_freq_eve(:,:,1,n_sample)')));
    end
    rate_k_para2 = rate_k_para2 / NumSamples_eve;
    
    rate_k_para3 = 0;
    for n_sample = 1:NumSamples
        rate_k_para3 = rate_k_para3 + log(abs(det(K(:,:,k) +  h_freq(:,:,1,n_sample,k) * diag(lambda(:,Nu+1)) * h_freq(:,:,1,n_sample,k)')));
    end
    rate_k_para3 = rate_k_para3 / NumSamples;
%    rate_k_para3 = log(abs(det(K(:,:,k) + cal_Gk_Aev_Gk(k,lambda,Omega))));
    rate_k_para4 = 0;
    for n_sample = 1:NumSamples_eve
        rate_k_para4 = rate_k_para4 + log(abs(det(eye(Ne) +  h_freq_eve(:,:,1,n_sample) * diag(lambda(:,k)+lambda(:,Nu+1)) * h_freq_eve(:,:,1,n_sample)')));
    end
    rate_k_para4 = rate_k_para4 / NumSamples_eve;
%    rate_k_para4 = log(abs(det(eye(Ne) +  cal_Gev_A_Gev(k,lambda,Omega,Omega_eve))));
%    rate_k_para4 = log(abs(det(eye(Ne) +  cal_Gev_A_Gev(k,lambda,Omega,Omega_eve))));
%    rate_all(k) = rate_k_para1 + rate_k_para2 - rate_k_para3 - rate_k_para4;
    rate_all(k) = rate_k_para1 + rate_k_para2 - rate_k_para3 - rate_k_para4; 
    rate_all(k) = max(rate_all(k),0);
end
sum_rate = sum(rate_all);
end

