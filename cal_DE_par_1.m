function [Gamma, Gamma_tilde,  Phi_tilde] = cal_DE_par_1(lambda, K_all, Omega)
% lambda  Nt * (Nu+1)
% K       Nr * Nr * Nu

[Nr,Nt,Nu] = size(Omega);
Phi_tilde = zeros(Nr,Nr,Nu);
Gamma = zeros(Nt,Nt,Nu);
Gamma_tilde = zeros(Nr,Nr,Nu);
for k = 1:Nu
    Phi_tilde_k = eye(Nr);
    while 1
        Phi_k = eye(Nt) + eta(Omega(:,:,k), inv(Phi_tilde_k)*inv(K_all(:,:,k))) * diag((lambda(:,k) + lambda(:,Nu+1)));
        Phi_tilde_k_new = eye(Nr) + eta_tilde(Omega(:,:,k), Phi_k\diag((lambda(:,k) + lambda(:,Nu+1)))) / K_all(:,:,k);
        error_diff_matrix = norm(Phi_tilde_k_new - Phi_tilde_k,'fro');
        if error_diff_matrix < 10^(-5)
            break;
        end
        Phi_tilde_k = Phi_tilde_k_new;
    end
    Phi_tilde(:,:,k) = Phi_tilde_k;
    Gamma(:,:,k) = eta(Omega(:,:,k),inv(Phi_tilde_k_new)*inv(K_all(:,:,k)));
    Gamma_tilde(:,:,k) = eta_tilde(Omega(:,:,k),Phi_k\diag((lambda(:,k) + lambda(:,Nu+1))));
end


end

