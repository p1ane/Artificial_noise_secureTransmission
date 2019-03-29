function [Out] = cal_Gev_A_Gev(k,lambda,Omega,Omega_ev)
%%%%%%%lambda  Nt * (Nu+1)
[Nrev,~] = size(Omega_ev);
[~,~,Nu] = size(Omega);

Out = zeros(Nrev,Nrev);

for n = 1:Nrev
    Out(n,n) = Omega_ev(n,:) * ( lambda(:,k) + lambda(:,Nu+1) );
end

end

