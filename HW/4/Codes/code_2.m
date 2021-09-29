%% 2 part 2
clc
clear
close all

signalSize = 1000;
M = 6;

%Sent signal 1:
mSource1 = 0;
sdSource1 = sqrt(10/2);
s1 = normrnd(mSource1, sdSource1, signalSize, 1) + 1i * normrnd(mSource1, sdSource1, signalSize, 1);

%Sent signal 2:
mSource2 = 0;
sdSource2 = sqrt(20/2);
s2 = normrnd(mSource2, sdSource2, signalSize, 1) + 1i * normrnd(mSource2, sdSource2, signalSize, 1);

%Noise info :
mNoise = 0;
sdNoise = sqrt(1/2);
noise = normrnd(mNoise, sdNoise, signalSize, M) + 1i * normrnd(mNoise, sdNoise, signalSize, M);

%For teta = 30: 
a30 = zeros(M,1);
for k=0:M-1
   a30(k+1) = exp(1i*pi*cos(pi/6)*k); 
end

%For teta = 160:
a160 = zeros(M,1);
for k=0:M-1
   a160(k+1) = exp(1i*pi*cos(8*pi/9)*k); 
end

%Recieved signal:
Ru = zeros(M,M);
for i=1:signalSize
    u = s1(i,1) * a30 + s2(i,1) * a160 + noise(i,:);
    Ru = Ru + u*u';
end
Ru = Ru / signalSize;

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

%% 2 part 3
clc
clear
close all

clc
clear
close all

signalSize = 100;
M = 6;

%Sent signal 1:
mSource1 = 0;
sdSource1 = sqrt(10/2);
s1 = normrnd(mSource1, sdSource1, signalSize, 1) + 1i * normrnd(mSource1, sdSource1, signalSize, 1);

%Noise info :
mNoise = 0;
sdNoise = sqrt(1/2);
noise = normrnd(mNoise, sdNoise, signalSize, M) + 1i * normrnd(mNoise, sdNoise, signalSize, M);

%For teta = 30: 
a30 = zeros(M,1);
for k=0:M-1
   a30(k+1) = exp(1i*pi*cos(pi/6)*k); 
end

%Input system:
u = s1 * a30.' + noise;
c = a30;
ca = null(a30.');

x = u * ca;
d = u * (c * inv(c'*c)); 

uPower = 10*sum(diag(c'*c))+1;
traceRx = 61*sum(diag(ca'*ca));
miu = 2 * 0.1 / traceRx;

wa = zeros(M-1,1);
for i=1:signalSize
    xvec = x(i,:);
    xvec = xvec.';
    e = d(i) - wa' * xvec;
    J(1,i) = e*conj(e);
    wa = wa + miu * e * xvec;
end

figure
n = 0:size(J,2)-1;
plot(n, J, 'Color', [147/255, 112/255, 219/255]);
