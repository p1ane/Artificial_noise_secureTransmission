function [Gra_UT_k] = cal_gra_to_Lambda(k,lambda,Omega,Omega_eve)
%%%%%%%%%%%%%%%%%%%%%%%%%lambda  Nt*(Nu+1)
[Nr,Nt,Nu] = size(Omega);
[Ne,~] = size(Omega_eve);

Gra_UT_k = zeros(Nt,Nu+1);    %%%%%%%%%%% Nu+1 是第k个用户减号项关于 Nu+1 个lambda求偏导
for lambda_i = 1:Nu+1           %%% UT_k关于lambda_i求导
    if lambda_i == k              %%% UT_k 关于 lambda_k 求导
        for m = 1:Nt
            temp = 0;
            for n = 1:Ne
                temp = temp + Omega_eve(n,m)/(1 + Omega_eve(n,:) * (lambda(:,k) + lambda(:,Nu+1)));
            end
            Gra_UT_k(m,lambda_i) = temp;
        end
        
    elseif lambda_i == Nu+1   %%% UT_k关于lambda_eve求导
        for m = 1:Nt
            temp_1 = 0;
            temp_2 = 0;
            for n = 1:Ne
                temp_1 = temp_1 + Omega_eve(n,m)/(1 + Omega_eve(n,:) * (lambda(:,k) + lambda(:,Nu+1)));
            end
            for n = 1:Nr
                lambda_cal = sum(lambda,2) - lambda(:,k);
                temp_2 = temp_2 + Omega(n,m,k)/(1 + Omega(n,:,k) * lambda_cal );
            end
            Gra_UT_k(m,lambda_i) = temp_1 + temp_2; 
        end
    else              %%lambda_i ~=k        %%% UT_k关于lambda求导 除了lambda_k 和 lambda_eve 
        for m = 1:Nt
            temp = 0;
            for n = 1:Nr
                lambda_cal = sum(lambda,2) - lambda(:,k);
                temp = temp + Omega(n,m,k)/(1 + Omega(n,:,k) * lambda_cal );
            end
            Gra_UT_k(m,lambda_i) = temp;
        end
    end

end

end

