%% 1 
%A

clc
clear 
close all 

P = [1;0.25];
R = [1.175 0.5; 0.5 1.175];
W = [0;0];
muOpt = 0.8511;
varD = 1;

[Wopt, J, Ws] = SDCalc(P, R, W, muOpt, varD, 2);

%with mu = 0.5 muOpt
[Wopt1, J1, Ws1] = SDCalc(P, R, W, 0.5 * muOpt, varD, 2);
%with mu = 0.1 muOpt
[Wopt2, J2, Ws2] = SDCalc(P, R, W, 0.1 * muOpt, varD, 2);

%with mu = muOpt
figure
n = 0:size(J,2)-1;
plot(n, J, 'Color', [147/255, 112/255, 219/255]);
title('Learning plot with optimum mu');
xlabel('n');
ylabel('J(w(n))')
hold on

%with mu = 0.5 muOpt
figure
n1 = 0:size(J1,2)-1;
plot(n1, J1, 'Color', [147/255, 112/255, 219/255]);
title('Learning plot with half optimum mu');
xlabel('n');
ylabel('J(w(n))')

%with mu = 0.1 muOpt
figure
n2 = 0:size(J2,2)-1;
plot(n2, J2, 'Color', [147/255, 112/255, 219/255]);
title('Learning plot with 0.1 optimum mu');
xlabel('n');
ylabel('J(w(n))')
hold off

%Comparing J with different mu 
figure
plot(n, J, 'Color', [147/255, 112/255, 219/255])
hold on 
plot(n1, J1, 'Color', [218/255, 112/255, 214/255]);
plot(n2, J2, 'Color', [102/255, 205/255, 170/255]);
title('Learning plot with different values of mu');
xlabel('n');
ylabel('J(w(n))')
legend('With optimum mu','With half optimum mu','With 0.1 optimum mu');
hold off


%%
%B - Plotting learning plot based on what been found and written in
%paper.

clc
clear
close all

%Finding Ru
Ru = zeros(9);
for i=1:9
    for j=1:9
        if i==j
            Ru(i,j) = 1.175;
        elseif abs(i-j) == 1
            Ru(i,j) = 0.5;
        elseif abs(i-j) == 2
            Ru(i,j) = 0.0625;
        end
    end
end

%Finding P
P = zeros(9,1);
P(4,1) = 0.25;
P(5,1) = 1;
P(6,1) = 0.25;

%Finding Wopt
Wopt = Ru\P;
disp('W optimum is : ');
disp(Wopt);

%Finding minimum squared error
Jmin = 1 - P' * inv(Ru) * P;

%Finding mu optimum
[eignVecRu, eigenValRu] = eig(Ru);
lambdaMax = max(eigenValRu(:));
lambdaMin = min(eigenValRu(eigenValRu>0));
muOpt = 2 / (lambdaMax + lambdaMin);

W = zeros(9,1);
varD = 1;

%Finding Wopt with simulation
[WoptSim, J, Ws] = SDCalc(P, Ru, W, muOpt, varD, 9);

%Plotting learning plot  
figure
n = 0:size(J,2)-1;
plot(n, J, 'Color', [147/255, 112/255, 219/255]);
title('Learning plot with optimum mu');
xlabel('n');
ylabel('J(w(n))')

%%
%2
clc
clear

clc
clear
close all

signalSize = 1000;
M = 5;

%Sent signal :
mSource = 0;
sdSource = sqrt(10);
s1 = normrnd(mSource, sdSource, signalSize, 1) + 1i * normrnd(mSource, sdSource, signalSize, 1);
s2 = normrnd(mSource, sdSource, signalSize, 1) + 1i * normrnd(mSource, sdSource, signalSize, 1);

%Noise info :
mNoise = 0;
sdNoise = 1;
noise = normrnd(mNoise, sdNoise, signalSize, 1) + 1i * normrnd(mNoise, sdNoise, signalSize, 1);

%For teta = 45 : 
a45 = zeros(5,1);
for k=0:M-1
   a45(k+1) = exp(1i*pi*cos(pi/4)*k); 
end

%For teta = 150
a150 = zeros(5,1);
for k=0:M-1
   a150(k+1) = exp(1i*pi*cos(5*pi/6)*k); 
end

%Recieved signal :
u = zeros(signalSize, M);
for j=1:M
    for i=1:signalSize
        u(i,j) = s1(i,1) * a45(j,1) + s2(i,1) * a150(j,1) + noise(i,1);
    end
end
u = u.';
%Findign Ru based on approximation given : 
Ru = RxApp(u, signalSize);

c45 = a45;
Wopt = (Ru\c45)*inv((c45')*(Ru\c45));

%%
%3
clc
clear
close all

signalSize = 1000;
M = 5;

%Sent signal :
mSource = 0;
sdSource = sqrt(5);
s1 = normrnd(mSource, sdSource, signalSize, 1) + 1i * normrnd(mSource, sdSource, signalSize, 1);
s2 = normrnd(mSource, sdSource, signalSize, 1) + 1i * normrnd(mSource, sdSource, signalSize, 1);

%Noise info :
mNoise = 0;
sdNoise = sqrt(0.5);
noise = normrnd(mNoise, sdNoise, signalSize, 1) + 1i * normrnd(mNoise, sdNoise, signalSize, 1);

%For teta = 45 : 
a45 = zeros(5,1);
for k=0:M-1
   a45(k+1) = exp(1i*pi*cos(pi/4)*k); 
end

%For teta = 150
a150 = zeros(5,1);
for k=0:M-1
   a150(k+1) = exp(1i*pi*cos(5*pi/6)*k); 
end

%Recieved signal :
u = zeros(signalSize, M);
for j=1:M
    for i=1:signalSize
        u(i,j) = s1(i,1) * a45(j,1) + s2(i,1) * a150(j,1) + noise(i,1);
    end
end
u = u.';
%Findign Ru based on approximation given : 
Ru = RxApp(u, signalSize);

%Defining teta and c(teta) :
if signalSize == 100
    teta = 0:(180/(signalSize)):180-1;
elseif signalSize == 1000
    teta = 0:(180/(signalSize+5)):180-1;
end    
teta = teta.';
cTeta = cTetaCal(teta, M);

%Jmin : 
Jmin = inv((cTeta')*(Ru\cTeta));
JminAvg = zeros(signalSize,1);
sum = sum(Jmin,2);
for i=1:signalSize
    JminAvg(i,1) = sum(i,1)/signalSize;
end
plot(abs(teta), abs(JminAvg), 'Color', [102/255, 205/255, 170/255])










