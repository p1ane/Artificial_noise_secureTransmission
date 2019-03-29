function [output] = eta(Omega,X)
output = diag(Omega'*diag(X));

end

