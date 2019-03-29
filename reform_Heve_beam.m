function [H] = reform_Heve_beam(Omega_eve,Nt,Nrev,NumSamples)
j = sqrt(-1);
H = zeros(Nrev,Nt,1,NumSamples);

    for sample_n = 1:NumSamples
        H(:,:,1,sample_n) = ...
        sqrt(0.5*Omega_eve(:,:)).*(randn(Nrev,Nt) + j*randn(Nrev,Nt));
    end



end

