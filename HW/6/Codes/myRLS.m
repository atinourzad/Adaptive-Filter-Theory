function [expJRLS] = myRLS(M, signalSize, iter, missAdj)
delta = 0.001;

lambda = (M - missAdj)/(missAdj + M);
fprintf('Calculated lambda is: %f \n ',lambda); 

J = zeros(signalSize, iter);
for a=1:iter
    
    %Making symbols:
    s = zeros(1,signalSize);
    symbols = [1/sqrt(2)+1i*1/sqrt(2), 1/sqrt(2)-1i*1/sqrt(2), -1/sqrt(2)+1i*1/sqrt(2), -1/sqrt(2)-1i*1/sqrt(2)];
    for i = 1:signalSize
        prob = rand(1);
        if prob < 0.25
            s(i) = symbols(1);
        elseif (prob > 0.25) &&  (prob < 0.5)
            s(i) = symbols(2);
        elseif (prob > 0.5) &&  (prob < 0.75)
            s(i) = symbols(3);
        else
            s(i) = symbols(4);
        end
    end

    %Creating u: (The output of channel)
    h = [0.5 1.2 1.5 -1];
    u = filter(h, 1, s);
 
    %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.01/2);
    noiseSize = signalSize;
    noise = normrnd(mNoise, sdNoise, 1, noiseSize) + 1i* normrnd(mNoise, sdNoise, 1, noiseSize);
    

    %Creating x: (A sum of channel output and noise)
    x = u + noise;

    %Creating the desired signal:
    d = zeros(1,signalSize);
    for i= 1:signalSize
        if i < 20
            d(i) = s(1);
        else
            d(i) = s(i-19);
        end
    end
    
    %RLS
    x = x.';
    w = zeros(M,1);
    p = eye(M)/delta;
    for n=M:signalSize
       xvec = x(n:-1:n-M+1);
       k = ((1/lambda) * p * xvec) / (1 + (1/lambda) * xvec' * p * xvec);
       alpha = d(n) - w' * xvec;
       w = w + k * conj(alpha);
       p = (1/lambda) * p - (1/lambda) * k * xvec' *  p;
       J(n,a) = alpha*conj(alpha);
    end
end
expJRLS = sum(J,2)/iter;
expJRLS(1:34)=[];
end

