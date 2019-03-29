function [Out] = cal_Gk_Aev_Gk(k,lambda,Omega)
%%%%%%%lambda  Nt * (Nu+1)
[Nr,~,Nu] = size(Omega);
Out = zeros(Nr,Nr);
for n = 1:Nr
    Out(n,n) = Omega(n,:,k) * lambda(:,Nu+1); 
end


end

