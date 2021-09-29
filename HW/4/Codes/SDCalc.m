function [J] = SDCalc(P, R, mu, varD, M, iter)
W = zeros(M,1);
J = 1;
for i=1:iter
    minuesGradJ = P - R * W;
    W = W + mu * minuesGradJ;
    J = JCalc(varD, P, W, R, J);
end
end

