function [Gamma, Gamma_tilde, Phi_tilde_new] = cal_DE_par_2(lambda,Omega,Omega_eve)
[Ne,Nt] = size(Omega_eve);
[~,~,Nu] = size(Omega);
% Gamma = zeros(Nt,Nt);
% Gamma_tilde = zeros(Ne,Ne);


Phi_tilde = eye(Ne);
while 1
    Phi = eye(Nt) + eta(Omega_eve, inv(Phi_tilde)) * diag(lambda(:,Nu+1));
    Phi_tilde_new = eye(Ne) + eta_tilde(Omega_eve, Phi\diag(lambda(:,Nu+1)));
    error_diff_matrix = norm(Phi_tilde_new - Phi_tilde,'fro');
    if error_diff_matrix < 10^(-5)
        break;
    end
    Phi_tilde = Phi_tilde_new;
end
Gamma = eta(Omega_eve, inv(Phi_tilde_new));
Gamma_tilde = eta_tilde(Omega_eve,Phi\diag(lambda(:,Nu+1)));
end

