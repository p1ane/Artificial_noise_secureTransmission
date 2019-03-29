function [K] = cal_K(lambda,Omega)
%%%%%%%%%%%%%%%%%%%%     I + \sum_{i\neq k}E{G_k*Q_i*G_K^H}
% lambda Nt * (Nu+1)     Nu个用户  1个窃听
% K      Nr * Nr * Nu
[Nr,Nt,Nu] = size(Omega);

for k = 1:Nu
    for antenna_i = 1:Nr
        temp_antenna = 0;
        for user_i = 1:Nu
            if user_i ~= k
                temp_antenna = temp_antenna + Omega(antenna_i,:,k) * lambda(:,user_i);
            end            
        end
        K_temp(antenna_i) = temp_antenna;
    end
    K(:,:,k) = diag(K_temp) + eye(Nr);
end


