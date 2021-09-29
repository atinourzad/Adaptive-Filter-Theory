%% 1 part 2
clc
clear
close all

M = 35;
signalSize = 4000;
iter = 100;

%Miu for M = 10%:
%The exact value :
% miu = 1.041168316e-3;
%The approximated value:
miu = 1.145285147e-3;


P = zeros(M, 1);
P(17,:) = -1;
P(18,:) = 1.5;
P(19,:) = 1.2;
P(20,:) = 0.5;

Rx = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Rx(i,j) = 4.9894;
        elseif abs(i-j) == 1
            Rx(i,j) = 0.9;
        elseif abs(i-j) == 2
            Rx(i,j) = -0.45;
        elseif abs(i-j) == 3
            Rx(i,j) = -0.5;
        end
    end
end

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
    u = zeros(1,signalSize);
    for k=1:signalSize
        u(1,k) = 0.5*s(1,k);
        if  k == 2
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1);
        elseif k == 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2);
        elseif k > 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2) - s(1,k-3);
        end
    end
 
    %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.0494/2);
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
expJ = sum(J,2)/iter;
expJ(1:34)=[];

Jmin = 1 - P' * inv(Rx) * P;
temp = 0;
for i=1:10
   temp = temp + expJ(size(expJ,1)-i);  
end
Jexcess = temp/10 - Jmin;
misadj = (Jexcess/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f',signalSize, iter, misadj); 

JSD = SDCalc(P, Rx, miu, 1, M, signalSize);
JminVec = ones(1,size(JSD,2))*Jmin;

%Plotting the learning figure
figure
nExpJ = 0:size(expJ,1)-1;
plot(nExpJ, expJ, 'Color', [147/255, 112/255, 219/255]);
hold on 
nJSD = 0:size(JSD,2)-1;
plot(nJSD, JSD, 'Color', [102/255, 205/255, 170/255]);
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);
title('Learning plot with optimum mu');
xlabel('n');
ylabel('E{J(w(n))}');
legend('LMS','SD','Jmin');
hold off

%% 1 part 3

clc
clear
close all

M = 35;
signalSize = 1000;
iter = 100;

%Miu for M = 10%:
%The exact value :
miu = 1.041168316e-3;
%The approximated value:
% miu = 1.145285147e-3;


P = zeros(M, 1);
P(17,:) = -1;
P(18,:) = 1.5;
P(19,:) = 1.2;
P(20,:) = 0.5;

Rx = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Rx(i,j) = 4.9894;
        elseif abs(i-j) == 1
            Rx(i,j) = 0.9;
        elseif abs(i-j) == 2
            Rx(i,j) = -0.45;
        elseif abs(i-j) == 3
            Rx(i,j) = -0.5;
        end
    end
end

for a=1:iter
    s = zeros(1,signalSize);
    symbols = [1/sqrt(2)+1i*1/sqrt(2), 1/sqrt(2)-1i*1/sqrt(2), -1/sqrt(2)+1i*1/sqrt(2), -1/sqrt(2)-1i*1/sqrt(2)];
    for i = 1:signalSize
        prob = rand(1);
        if prob < 0.25 || prob == 0.25
            s(i) = symbols(1);
        elseif (prob > 0.25) &&  (prob < 0.5)
            s(i) = symbols(2);
        elseif (prob > 0.5) &&  (prob < 0.75)
            disp('hi');
            s(i) = symbols(3);
        else
            s(i) = symbols(4);
        end
    end

    %Creating u: (The output of channel)
    u = zeros(1,signalSize);
    for k=1:signalSize
        u(1,k) = 0.5*s(1,k);
        if  k == 2
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1);
        elseif k == 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2);
        elseif k > 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2) - s(1,k-3);
        end
    end
 
    %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.0494/2);
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
    
    %LMS
    x = x.';
    w = zeros(M,1);
    y = zeros(iter,signalSize-M+1);
    for n=M:signalSize
       xvec = x(n:-1:n-M+1);
       y(a,n-M+1) = w' * xvec;
       e = d(n) - w' * xvec;
       w = w + miu * xvec * conj(e);
    end
end
yNew = sum(y,1);

figure
scatter(real(x), imag(x), 'filled');
figure
c = linspace(1,5,length(real(yNew)));
scatter(real(yNew), imag(yNew),[], c, 'filled');
%% 1 part 4 M = 5%
clc
clear
close all

M = 35;
signalSize = 4000;
iter = 100;

%Miu for M = 5%:
%The exact value :
% miu = 5.453739e-4;
%The approximated value:
miu = 5.726425e-4;

P = zeros(M, 1);
P(17,:) = -1;
P(18,:) = 1.5;
P(19,:) = 1.2;
P(20,:) = 0.5;

Rx = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Rx(i,j) = 4.9894;
        elseif abs(i-j) == 1
            Rx(i,j) = 0.9;
        elseif abs(i-j) == 2
            Rx(i,j) = -0.45;
        elseif abs(i-j) == 3
            Rx(i,j) = -0.5;
        end
    end
end

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
    u = zeros(1,signalSize);
    for k=1:signalSize
        u(1,k) = 0.5*s(1,k);
        if  k == 2
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1);
        elseif k == 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2);
        elseif k > 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2) - s(1,k-3);
        end
    end
 
    %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.0494/2);
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
expJ = sum(J,2)/iter;
expJ(1:34)=[];

Jmin = 1 - P' * inv(Rx) * P;
temp = 0;
for i=1:10
   temp = temp + expJ(size(expJ,1)-i);  
end
Jexcess = temp/10 - Jmin;
misadj = (Jexcess/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f',signalSize, iter, misadj); 

JSD = SDCalc(P, Rx, miu, 1, M, signalSize);
JminVec = ones(1,size(JSD,2))*Jmin;

%Plotting the learning figure
figure
nExpJ = 0:size(expJ,1)-1;
plot(nExpJ, expJ, 'Color', [147/255, 112/255, 219/255]);
hold on 
nJSD = 0:size(JSD,2)-1;
plot(nJSD, JSD, 'Color', [102/255, 205/255, 170/255]);
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);
title('Learning plot with optimum mu');
xlabel('n');
ylabel('E{J(w(n))}');
legend('LMS','SD','Jmin');
hold off

%% 1 part 4 M = 1%
clc
clear
close all

M = 35;
signalSize = 30000;
iter = 100;

%Miu for M = 1%:
%The exact value :
% miu = 1.1339e-4;
%The approximated value:
miu = 1.145285e-4;


P = zeros(M, 1);
P(17,:) = -1;
P(18,:) = 1.5;
P(19,:) = 1.2;
P(20,:) = 0.5;

Rx = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Rx(i,j) = 4.9894;
        elseif abs(i-j) == 1
            Rx(i,j) = 0.9;
        elseif abs(i-j) == 2
            Rx(i,j) = -0.45;
        elseif abs(i-j) == 3
            Rx(i,j) = -0.5;
        end
    end
end

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
    u = zeros(1,signalSize);
    for k=1:signalSize
        u(1,k) = 0.5*s(1,k);
        if  k == 2
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1);
        elseif k == 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2);
        elseif k > 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2) - s(1,k-3);
        end
    end
 
    %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.0494/2);
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
expJ = sum(J,2)/iter;
expJ(1:34)=[];

Jmin = 1 - P' * inv(Rx) * P;
temp = 0;
for i=1:10
   temp = temp + expJ(size(expJ,1)-i);  
end
Jexcess = temp/10 - Jmin;
misadj = (Jexcess/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f',signalSize, iter, misadj); 

JSD = SDCalc(P, Rx, miu, 1, M, signalSize);
JminVec = ones(1,size(JSD,2))*Jmin;

%Plotting the learning figure
figure
nExpJ = 0:size(expJ,1)-1;
plot(nExpJ, expJ, 'Color', [147/255, 112/255, 219/255]);
hold on 
nJSD = 0:size(JSD,2)-1;
plot(nJSD, JSD, 'Color', [102/255, 205/255, 170/255]);
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);
title('Learning plot with optimum mu');
xlabel('n');
ylabel('E{J(w(n))}');
legend('LMS','SD','Jmin');
hold off

%% 1 part 5
clc
clear
close all

M = 35;
signalSize = 4000;
iter = 100;
symbols = [1/sqrt(2)+1i*1/sqrt(2), 1/sqrt(2)-1i*1/sqrt(2), -1/sqrt(2)+1i*1/sqrt(2), -1/sqrt(2)-1i*1/sqrt(2)];
%Miu for M = 10%:
%The exact value :
% miu = 1.041168316e-3;
%The approximated value:
miu = 1.145285147e-3;


P = zeros(M, 1);
P(17,:) = -1;
P(18,:) = 1.5;
P(19,:) = 1.2;
P(20,:) = 0.5;

Rx = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Rx(i,j) = 4.9894;
        elseif abs(i-j) == 1
            Rx(i,j) = 0.9;
        elseif abs(i-j) == 2
            Rx(i,j) = -0.45;
        elseif abs(i-j) == 3
            Rx(i,j) = -0.5;
        end
    end
end

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
    u = zeros(1,signalSize);
    for k=1:signalSize
        u(1,k) = 0.5*s(1,k);
        if  k == 2
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1);
        elseif k == 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2);
        elseif k > 3
            u(1,k) = 0.5*s(1,k) + 1.2 * s(1,k-1) + 1.5 * s(1,k-2) - s(1,k-3);
        end
    end
 
    %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.0494/2);
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
    
    %LMS
    x = x.';
    w = zeros(M,1);
    dLearnt = zeros(1,signalSize);
    for n=M:signalSize
        if n < 55 
            xvec = x(n:-1:n-M+1);
            y = w' * xvec;
            if real(y) > 0 && imag(y) > 0 
                dLearnt(1,n) = symbols(1);
            elseif real(y) > 0 && imag(y) < 0
                dLearnt(1,n) = symbols(2);
            elseif real(y) < 0 && imag(y) > 0
                dLearnt(1,n) = symbols(3);
            else
                dLearnt(1,n) = symbols(4);
            end
            e = d(n) - w' * xvec;
            J(n,a) = e*conj(e);
            w = w + miu * xvec * conj(e);
        else
            xvec = x(n:-1:n-M+1);
            y = w' * xvec;
            if real(y) > 0 && imag(y) > 0 
                dLearnt(1,n) = symbols(1);
            elseif real(y) > 0 && imag(y) < 0
                dLearnt(1,n) = symbols(2);
            elseif real(y) < 0 && imag(y) > 0
                dLearnt(1,n) = symbols(3);
            else
                dLearnt(1,n) = symbols(4);
            end
            e = dLearnt(1,n) - w' * xvec;
            J(n,a) = e*conj(e);
            w = w + miu * xvec * conj(e);
        end
    end
end

expJ = sum(J,2)/iter;
expJ(1:34)=[];

Jmin = 1 - P' * inv(Rx) * P;
temp = 0;
for i=1:10
   temp = temp + expJ(size(expJ,1)-i);  
end
Jexcess = temp/10 - Jmin;
misadj = (Jexcess/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f',signalSize, iter, misadj); 

JSD = SDCalc(P, Rx, miu, 1, M, signalSize);
JminVec = ones(1,size(JSD,2))*Jmin;

%Plotting the learning figure
figure
nExpJ = 0:size(expJ,1)-1;
plot(nExpJ, expJ, 'Color', [147/255, 112/255, 219/255]);
hold on 
nJSD = 0:size(JSD,2)-1;
plot(nJSD, JSD, 'Color', [102/255, 205/255, 170/255]);
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);
title('Learning plot with optimum mu');
xlabel('n');
ylabel('E{J(w(n))}');
legend('LMS','SD','Jmin');
hold off

        