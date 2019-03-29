function [object] = cvx_object(lambda_new,lambda_i,Omega,Omega_eve)
[Nr,Nt,Nu] = size(Omega);
[Ne,~] = size(Omega_eve);
[K_i] = cal_K(lambda_i,Omega);
[K_new] = cal_K(lambda_new,Omega);
[Gamma_1,Gamma_tilde_1,Phi_tilde_1] = cal_DE_par_1(lambda_i,K_i,Omega);
[Gamma_2,Gamma_tilde_2,Phi_tilde_2] = cal_DE_par_2(lambda_i,Omega,Omega_eve);
object1 = 0;      %%%%%%%% 第一项加法项 A + B
object2 = 0;      %%%%%%%% 第二项减法项 C + D
DE_2_temp = log_det(eye(Nt) + Gamma_2 * diag(lambda_new(:,Nu+1)))...
                      + log_det(Gamma_tilde_2 + eye(Ne));    %%%%% 确定性等同第二项，每个用户都一样，只算一次
for k = 1:Nu
    object1 = object1 + log_det(eye(Nt)+Gamma_1(:,:,k) * diag(lambda_new(:,k) + lambda_new(:,Nu+1)))...
                      + log_det(Gamma_tilde_1(:,:,k) + K_new(:,:,k) );    
end

object1 = object1 +  Nu * DE_2_temp;

for k = 1:Nu
    [Gra_UT_k] = cal_gra_to_Lambda(k,lambda_i,Omega,Omega_eve);
    for kk = 1:Nu + 1
      object2 = object2 + Gra_UT_k(:,kk)' * (lambda_new(:,kk) - lambda_i(:,kk)); 
    end
end

object = object1 - object2; 



end

