%% Part 1
clc
clear
close all

M = 35;
signalSize = 4000;
iter = 100;

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

missAdjTheory = 0.1;
miu = 2*missAdjTheory/trace(Rx);
fprintf('For given miss adjustment the calculated miu is :%f \n ',miu); 


for a=1:iter
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
    
    %LMS
    x = x.';
    w = zeros(M,1);
    for n=M:signalSize
       xvec = x(n:-1:n-M+1);
       e = d(n) - w' * xvec;
       J(n,a) = e*conj(e);
       w = w + miu * xvec * conj(e);
    end
end
expJLMS = sum(J,2)/iter;
expJLMS(1:34)=[];

Jmin = 1 - P' * inv(Rx) * P;
temp = 0;
for i=1:100
   temp = temp + expJLMS(size(expJLMS,1)-i);  
end
Jexcess = temp/100 - Jmin;
misadj = (Jexcess/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f \n ',signalSize, iter, misadj); 

JminVec = ones(1,size(expJLMS,1))*Jmin;

%Plotting the learning curve
figure
nExpJ = 0:size(expJLMS,1)-1;
plot(nExpJ, expJLMS, 'Color', [147/255, 112/255, 219/255]);
hold on 
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);
title('Learning curve with optimum mu');
xlabel('n');
ylabel('E{J(w(n))}');
legend('LMS','Jmin');
hold off

