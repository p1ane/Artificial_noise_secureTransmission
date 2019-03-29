function [output] = eta_tilde(Omega,X)
output = diag(Omega*diag(X));
end

