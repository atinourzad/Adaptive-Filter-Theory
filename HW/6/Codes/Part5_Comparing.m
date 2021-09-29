%% Part 5 Comparing

clc
close all
clear

M = 35;
signalSize = 700;
iter = 100;
delta = 798;

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

Jmin = 1 - P' * inv(Rx) * P;

missAdjTheory = 10/100;
lambda = (M - missAdjTheory)/(missAdjTheory + M);
fprintf('Calculated lambda for M =%d and miss adjustment 0.1 is : %f \n ',M ,lambda); 

miu = 2*missAdjTheory/trace(Rx);
fprintf('For given miss adjustment the calculated miu is :%f \n ',miu);

JavgRLSL = zeros(M,signalSize);
JLMS = zeros(signalSize,iter);
JRLS = zeros(signalSize,iter);
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
    x = x.';
    
    %LMS
    wLMS = zeros(M,1);
    for n=M:signalSize
       xvecLMS = x(n:-1:n-M+1);
       eLMS = d(n) - wLMS' * xvecLMS;
       JLMS(n,a) = eLMS*conj(eLMS);
       wLMS = wLMS + miu * xvecLMS * conj(eLMS);
    end
    
    %RLS
    wRLS = zeros(M,1);
    deltaRLS = 0.001;
    pRLS = eye(M)/deltaRLS;
    for n=M:signalSize
       xvec = x(n:-1:n-M+1);
       k = ((1/lambda) * pRLS * xvec) / (1 + (1/lambda) * xvec' * pRLS * xvec);
       alpha = d(n) - wRLS' * xvec;
       wRLS = wRLS + k * conj(alpha);
       pRLS = (1/lambda) * pRLS - (1/lambda) * k * xvec' *  pRLS;
       JRLS(n,a) = alpha*conj(alpha);
    end
    
    %RLSL
    gamma = ones(M+1, signalSize);

    Jf = (delta)*(lambda^(-2))*ones(M,signalSize);
    Jb = ones(M,signalSize);
    for i=1:signalSize
        Jb(:,i) = (delta) * (lambda.^(-2:-1:-M-1));
    end
    Jfb = zeros(M,signalSize);
    Jeb = zeros(M,signalSize);

    b = zeros(M+1,signalSize);
    f = zeros(M+1,signalSize);
    eLMS = zeros(M+1,signalSize);

    kf = zeros(M,signalSize);
    kb = zeros(M,signalSize);
    h = zeros(M,signalSize);

    for n=2:signalSize
        gamma(1,n-1) = 1;
        f(1,n) = x(n);
        b(1,n) = x(n);
        eLMS(1,n) = d(n);

        for m = 1:M
            Jf(m,n) = lambda * Jf(m,n-1) + abs(f(m,n))^2 / gamma(m,n-1);
            Jb(m,n) = lambda * Jb(m,n-1) + abs(b(m,n))^2 / gamma(m,n);

            Jfb(m,n) = lambda * Jfb(m,n-1) + (conj(b(m,n-1))*f(m,n)) / gamma(m,n-1);
            Jeb(m,n) = lambda * Jeb(m,n-1) + (eLMS(m,n)*conj(b(m,n))) / gamma(m,n);

            kf(m+1,n) = -conj(Jfb(m,n)) / Jb(m,n-1);
            kb(m+1,n) = -Jfb(m,n) / Jf(m,n);
            h(m,n) = conj(Jeb(m,n))/Jb(m,n);

            f(m+1,n) = f(m,n) + conj(kf(m+1,n)) * b(m,n-1);
            b(m+1,n) = b(m,n-1) + conj(kb(m+1,n)) * f(m,n);
            eLMS(m+1,n) = eLMS(m,n) - conj(h(m,n)) * b(m,n);

            gamma(m+1,n) = gamma(m,n) - ((abs(b(m,n))^2)/(Jb(m,n)));
        end 
    end
    JavgRLSL = JavgRLSL + abs(eLMS(1:end-1,:)).^2;
end
JavgRLSL = JavgRLSL/iter;

expJLMS = sum(JLMS,2)/iter;
expJLMS(1:34)=[];

expJRLS = sum(JRLS,2)/iter;
expJRLS(1:34)=[];

%Plotting the learning curve
figure
JavgNewTap35 = JavgRLSL(35,:);
nExpJ35 = 0:size(JavgNewTap35,2)-1;
plot(nExpJ35, JavgNewTap35, 'Color', [102/255, 205/255, 170/255]);

hold on 
nExpJ = 0:size(expJRLS,1)-1;
plot(nExpJ, expJRLS, 'Color', [255/255, 0/255, 255/255]);

nExpJ = 0:size(expJLMS,1)-1;
plot(nExpJ, expJLMS, 'Color', [147/255, 112/255, 219/255]);

JminVec = ones(1,size(expJLMS,1))*Jmin;
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);

title('Learning curve with optimum mu');
xlabel('n');
ylabel('E{J}');
legend('RLSL','RLS', 'LMS', 'Jmin');
hold off

