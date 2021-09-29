%% Part 5

clc
clear 
close all

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
fprintf('Calculated lambda for M =%d and miss adjustment %10 is : %f \n ',M ,lambda); 

M = 50;
Javg = zeros(M,signalSize);
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
    
    %RLSL
    x = x.';
    
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
    e = zeros(M+1,signalSize);

    kf = zeros(M,signalSize);
    kb = zeros(M,signalSize);
    h = zeros(M,signalSize);

    for n=2:signalSize
        gamma(1,n-1) = 1;
        f(1,n) = x(n);
        b(1,n) = x(n);
        e(1,n) = d(n);

        for m = 1:M
            Jf(m,n) = lambda * Jf(m,n-1) + abs(f(m,n))^2 / gamma(m,n-1);
            Jb(m,n) = lambda * Jb(m,n-1) + abs(b(m,n))^2 / gamma(m,n);

            Jfb(m,n) = lambda * Jfb(m,n-1) + (conj(b(m,n-1))*f(m,n)) / gamma(m,n-1);
            Jeb(m,n) = lambda * Jeb(m,n-1) + (e(m,n)*conj(b(m,n))) / gamma(m,n);

            kf(m+1,n) = -conj(Jfb(m,n)) / Jb(m,n-1);
            kb(m+1,n) = -Jfb(m,n) / Jf(m,n);
            h(m,n) = conj(Jeb(m,n))/Jb(m,n);

            f(m+1,n) = f(m,n) + conj(kf(m+1,n)) * b(m,n-1);
            b(m+1,n) = b(m,n-1) + conj(kb(m+1,n)) * f(m,n);
            e(m+1,n) = e(m,n) - conj(h(m,n)) * b(m,n);

            gamma(m+1,n) = gamma(m,n) - ((abs(b(m,n))^2)/(Jb(m,n)));
        end 
    end
    Javg = Javg + abs(e(1:end-1,:)).^2;
end
Javg = Javg/iter;

JavgNewTap = Javg(35,:);
Jexcess = sum(JavgNewTap(end-100:end))/100-Jmin;
missAdj = Jexcess * 100 /Jmin;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, miss adjustment is :%f \n ',signalSize, iter, missAdj); 

%Plotting the learning curve
figure

JavgNewTap30 = Javg(30,:);
nExpJ30 = 0:size(JavgNewTap30,2)-1;
plot(nExpJ30, JavgNewTap30, 'Color', [147/255, 112/255, 219/255]);

hold on 
JavgNewTap35 = Javg(35,:);
nExpJ35 = 0:size(JavgNewTap35,2)-1;
plot(nExpJ35, JavgNewTap35, 'Color', [102/255, 205/255, 170/255]);

JavgNewTap40 = Javg(40,:);
nExpJ40 = 0:size(JavgNewTap40,2)-1;
plot(nExpJ40, JavgNewTap40, 'Color', [255/255, 0/255, 255/255]);

JavgNewTap45 = Javg(45,:);
nExpJ45 = 0:size(JavgNewTap45,2)-1;
plot(nExpJ45, JavgNewTap45, 'Color', [255/255, 215/255, 0/255]);

JavgNewTap50 = Javg(50,:);
nExpJ50 = 0:size(JavgNewTap50,2)-1;
plot(nExpJ50, JavgNewTap50, 'Color', [0/255, 191/255, 255/255]);

JminVec = ones(1,signalSize)*Jmin;
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [0/255, 0/255, 255/255]);

title('Learning curve of RLSL');
xlabel('n');
ylabel('E{J(w(n))}');
legend('M=30','M=35','M=40','M=45', 'M=50', 'Jmin');
hold off

