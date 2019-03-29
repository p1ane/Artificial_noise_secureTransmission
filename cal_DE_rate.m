function [sum_rate, rate_all] = cal_DE_rate(lambda,Omega,Omega_eve)
[Nr,Nt,Nu] = size(Omega);
[Ne,~] = size(Omega_eve);
[K] = cal_K(lambda,Omega);

rate_all = zeros(Nu,1);
[Gamma_1, Gamma_tilde_1,  Phi_tilde_1] = cal_DE_par_1(lambda, K, Omega);
[Gamma_2, Gamma_tilde_2,  Phi_tilde_2] = cal_DE_par_2(lambda, Omega, Omega_eve);
for k = 1:Nu
    
    rate_k_para1 = log(abs(det(eye(Nt) + Gamma_1(:,:,k) * diag( (lambda(:,k)+lambda(:,Nu+1))))))...
        + log(abs(det(Gamma_tilde_1(:,:,k)+K(:,:,k)))) - trace(eye(Nr) - inv(Phi_tilde_1(:,:,k)));
    rate_k_para2 = log(abs(det(eye(Nt) +  Gamma_2 * diag(lambda(:,Nu+1)) ))) ...
        + log(abs(det(Gamma_tilde_2 + eye(Ne)))) - trace(eye(Ne) - inv(Phi_tilde_2));

    rate_k_para3 = log(abs(det(K(:,:,k) + cal_Gk_Aev_Gk(k,lambda,Omega))));
    rate_k_para4 = log(abs(det(eye(Ne) +  cal_Gev_A_Gev(k,lambda,Omega,Omega_eve))));
    
    rate_all(k) = rate_k_para1 + rate_k_para2 - rate_k_para3 - rate_k_para4;
    
    rate_all(k) = max(rate_all(k),0);
end
sum_rate = sum(rate_all);

end

