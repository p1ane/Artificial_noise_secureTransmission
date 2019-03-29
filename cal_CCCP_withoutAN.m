function [lambda_optimal,rate_iter] = cal_CCCP_withoutAN(lambda,P,Omega,Omega_eve,Looptime)
[Nr,Nt,Nu] = size(Omega);
[Ne,~] = size(Omega_eve);
lambda_i = lambda;
rate_i = 0;
rate_iter = zeros(1,Looptime);
for loop = 1:Looptime
    cvx_begin sdp
    variable lambda_new(Nt,Nu+1)
    expression object
    object = cvx_object(lambda_new,lambda_i,Omega,Omega_eve);
    maximize object
    subject to
    sum(sum(lambda_new)) <= P
    lambda_new(:,Nu+1) == 0;
    for i = 1:Nu+1
        0 <= lambda_new(:,i);
    end

    cvx_end
    
    [rate_new,~] = cal_DE_rate(lambda_new,Omega,Omega_eve);
%     if abs(rate_new - rate_i) / rate_i < 5e-4
%         break;
%     end
    lambda_i = lambda_new;
    rate_i = rate_new
    rate_iter(loop) = rate_new;
    
end

lambda_optimal = lambda_new;

end

