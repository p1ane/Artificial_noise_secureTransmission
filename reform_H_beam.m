function [H] = reform_H_beam(Omega,Nt,Nr,NumSamples,Nu)
j = sqrt(-1);
H = zeros(Nr,Nt,1,NumSamples,Nu);
for k = 1:Nu
    for sample_n = 1:NumSamples
        H(:,:,1,sample_n,k) = ...
        sqrt(0.5*Omega(:,:,k)).*(randn(Nr,Nt) + j*randn(Nr,Nt));
    end
end


end

