function [W, J, Ws] = SDCalc(P, R, W, mu, varD, M)
n=1;
endOfLoop = zeros(M,1);
error = 1e-5;
appEq = @(x,y,err) abs(x-y)<err;
Ws = [0;0];
J = 1;
while 1
    minuesGradJ = P - R * W;
    endFlag = appEq(minuesGradJ, endOfLoop, error);
 
    if endFlag
        break
    else
        W = W + mu * minuesGradJ;
        J = JCalc(varD, P, W, R, J);
        n = n + 1;
        Ws(1,n) = W(1,1);
        Ws(2,n) = W(2,1);
    end
end
end

