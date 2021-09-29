%% Part 2

clc
clear
close all

M = 35;
signalSize = 1000;
iter = 100;
delta = 0.001;
missAdjTheory = 0.1;

lambda = (M - missAdjTheory)/(missAdjTheory + M);
fprintf('Calculated lambda is: %f \n ',lambda); 

P = zeros(M, 1);
P(17,:) = -1;
P(18,:) = 1.5;
P(19,:) = 1.2;
P(20,:) = 0.5;

Rx = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Rx(i,j) = 4.95;
        elseif abs(i-j) == 1
            Rx(i,j) = 0.9;
        elseif abs(i-j) == 2
            Rx(i,j) = -0.45;
        elseif abs(i-j) == 3
            Rx(i,j) = -0.5;
        end
    end
end

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
        if i < 18
            d(i) = s(1);
        else
            d(i) = s(i-17);
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

Jmin = 1 - P' * inv(Rx) * P;

temp = 0;
for i=1:100
   temp = temp + expJRLS(size(expJRLS,1)-i);  
end

Jexcess = temp/100 - Jmin;
misadj = (Jexcess/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f \n ',signalSize, iter, misadj); 

JminVec = ones(1,size(expJRLS,1))*Jmin;

%Plotting the learning curve
figure
nExpJ = 0:size(expJRLS,1)-1;
plot(nExpJ, expJRLS, 'Color', [147/255, 112/255, 219/255]);
hold on 
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);
title('Learning curve with optimum mu');
xlabel('n');
ylabel('E{J}');
legend('RLS','Jmin');
hold off
