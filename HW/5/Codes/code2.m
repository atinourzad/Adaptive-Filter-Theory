%% 1 part 1
clc
clear
close all

M = 8;
signalSize = 4000;
iter = 100;

W = [0.4 1 -0.3];
H = [0.1 -0.2 0.5 1 -1 -0.5 0.3];

Ru = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Ru(i,j) = 10.56;
        elseif abs(i-j) == 1
            Ru(i,j) = -1.08;
        elseif abs(i-j) == 2
            Ru(i,j) = -5.8;
        elseif abs(i-j) == 3
            Ru(i,j) = 1.4;
        elseif abs(i-j) == 4
            Ru(i,j) = 0.6;
        elseif abs(i-j) == 5
            Ru(i,j) = -0.44;
        elseif abs(i-j) == 6
            Ru(i,j) = 0.12;
        elseif abs(i-j) == 7
            Ru(i,j) = 0;
        end
    end
end

missAdjTheory = 0.05;
miu = 2*missAdjTheory/trace(Ru);

P = zeros(M, 1);
for i=1:M
    if i == 1
        P(i,1) = 0.4 * Ru(1,i) + Ru(1,i+1) - 0.3 * Ru(1,i+2);
    elseif i == 2
        P(i,1) = 0.4 * Ru(1,i) + Ru(1,i-1) - 0.3 * Ru(1,i);
    elseif i > 2
        P(i,1) = 0.4 * Ru(1,i) + Ru(i-1,1) - 0.3 * Ru(1,i-2);
    end
end

for i = 1:iter
    
   %Creating u[n]
   mSignal = 0;
   sdSignal = sqrt(4);
   input = normrnd(mSignal, sdSignal, 1, signalSize);
   
   u = filter(H, 1, input);
   
   %Creating uNew : output of w filter
   uNew = filter(W, 1, u);
   
   %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.05);
    noiseSize = signalSize;
    noise = normrnd(mNoise, sdNoise, 1, noiseSize); 
    
    %Creating x: result of adding white noise to uNew
    x = uNew + noise;
    
    %Performing LMS:
    u = u.';
    x = x.';
    w = zeros(M,1);
    for n=M:signalSize
       uvec = u(n:-1:n-M+1);
       e = x(n) - w' * uvec;
       J(n,i) = e*conj(e);
       w = w + miu * uvec * conj(e);
    end
end

expJLMS = sum(J,2)/iter;
expJLMS(1:M-1)=[];

Jmin = 14.426 - P' * inv(Ru) * P;

temp = 0;
for i=1:100
   temp = temp + expJLMS(size(expJLMS,1)-i);  
end
Jexcess = temp/100 - Jmin;
misadj = (Jexcess/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f \n ',signalSize, iter, misadj); 

JminVec = ones(1,size(expJLMS,1))*Jmin;

%Plotting the learning figure
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

%Comparing w of LMS with W:
figure
freqz(W,1, 'whole');
hold on
freqz(w,1, 'whole');
title('Comparing adaptive filter weigths with the actual weights');
legend('Actual filter Weights','Adaptive Filter Weights');
lines = findall(gcf,'type','line');
lines(1).Color = [147/255, 112/255, 219/255];
lines(2).Color = [255/255, 215/255, 0/255];
lines(3).Color = [147/255, 112/255, 219/255];
hold off

% finding spread value:
eignValues = eig(Ru);
lambdaMax = max(eignValues);
lambdaMin = min(eignValues);
spreadValue = lambdaMin / lambdaMax;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, spread value is :%f \n ',signalSize, iter, spreadValue); 
%% 2 part 2

clc
clear
close all

M = 8;
signalSize = 4000;
iter = 100;
beta = 0.9;

W = [0.4 1 -0.3];
H = [0.1 -0.2 0.5 1 -1 -0.5 0.3];

Ru = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Ru(i,j) = 10.56;
        elseif abs(i-j) == 1
            Ru(i,j) = -1.08;
        elseif abs(i-j) == 2
            Ru(i,j) = -5.8;
        elseif abs(i-j) == 3
            Ru(i,j) = 1.4;
        elseif abs(i-j) == 4
            Ru(i,j) = 0.6;
        elseif abs(i-j) == 5
            Ru(i,j) = -0.44;
        elseif abs(i-j) == 6
            Ru(i,j) = 0.12;
        elseif abs(i-j) == 7
            Ru(i,j) = 0;
        end
    end
end

missAdjTheory = 0.05;
miu = 2*missAdjTheory/trace(Ru);
miuHat = 2*missAdjTheory/M;

P = zeros(M, 1);
for i=1:M
    if i == 1
        P(i,1) = 0.4 * Ru(1,i) + Ru(1,i+1) - 0.3 * Ru(1,i+2);
    elseif i == 2
        P(i,1) = 0.4 * Ru(1,i) + Ru(1,i-1) - 0.3 * Ru(1,i);
    elseif i > 2
        P(i,1) = 0.4 * Ru(1,i) + Ru(i-1,1) - 0.3 * Ru(1,i-2);
    end
end

%DCT Transform:
TDCT = zeros(M,M);
for i=1:M
   if i == 1
       TDCT(1,:) = 1/sqrt(M);
   else
       for j=1:M
          TDCT(i,j) = sqrt(2/M)*cos((pi*(2*j-1)*(i-1))/(2*M)); 
       end
   end
end

for i = 1:iter
    
   %Creating u[n]
   mSignal = 0;
   sdSignal = sqrt(4);
   input = normrnd(mSignal, sdSignal, 1, signalSize);
   u = filter(H, 1, input);
   
   %Creating uNew : output of w filter
   uNew = filter(W, 1, u);
   
   %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.05);
    noiseSize = signalSize;
    noise = normrnd(mNoise, sdNoise, 1, noiseSize); 
    
    %Creating x: result of adding white noise to uNew
    x = uNew + noise;
    
    %Performing LMS:
    u = u.';
    x = x.';
    wLMS = zeros(M,1);
    for n=M:signalSize
       uvec = u(n:-1:n-M+1);
       e = x(n) - wLMS' * uvec;
       JLMS(n,i) = e*conj(e);
       wLMS = wLMS + miu * uvec * conj(e);
    end
    
    %Performing DCT
    wDCT = zeros(M,1);
    sigmaDCT = ones(M,M);
    for n=M:signalSize
       uvec = u(n:-1:n-M+1);
       uVec = TDCT * uvec;
       e = x(n) - wDCT' * uVec;
       for k=1:M
           wDCT(k) = wDCT(k) + miuHat * uVec(k) * conj(e) / sigmaDCT(k,k);
           sigmaDCT(k,k) = beta * sigmaDCT(k,k) + (1-beta) * (abs(uVec(k)).^2);
       end
       JDCT(n,i) = e*conj(e);
    end
end

expJLMS = sum(JLMS,2)/iter;
expJLMS(1:M-1)=[];

expJDST = sum(JDCT,2)/iter;
expJDST(1:M-1)=[];

Jmin = 14.426 - P' * inv(Ru) * P;
JminVec = ones(1,size(expJLMS,1))*Jmin;


temp = 0;
for i=1:100
   temp = temp + expJLMS(size(expJLMS,1)-i);  
end
JexcessLMS = temp/100 - Jmin;
misAdjLMS = (JexcessLMS/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment of LMS is :%f \n ',signalSize, iter, misAdjLMS); 
newline();
temp = 0;
for i=1:100
   temp = temp + expJDST(size(expJDST,1)-i);  
end
JexcessDFT = temp/100 - Jmin;
misAdjDCT = (JexcessDFT/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment of DCT is :%f \n ',signalSize, iter, misAdjDCT); 
newline();
%Plotting the learning figure
figure
nExpJLMS = 0:size(expJLMS,1)-1;
plot(nExpJLMS, expJLMS, 'Color', [147/255, 112/255, 219/255]);
hold on 
nExpJDCT = 0:size(expJDST,1)-1;
plot(nExpJDCT, expJDST, 'Color', [102/255, 205/255, 170/255]);
nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);

title('Learning curve with optimum mu');
xlabel('n');
ylabel('E{J(w(n))}');
legend('LMS','DCT','Jmin');
hold off

% finding spread value:
RuVec = TDCT * Ru * TDCT;
eignValues = eig(sigmaDCT^(-0.5)*RuVec*sigmaDCT^(-0.5));
lambdaMax = max(eignValues);
lambdaMin = min(eignValues);
spreadValue = lambdaMin / lambdaMax;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, spread value is :%f \n ',signalSize, iter, spreadValue); 
%% 2 part 3

clc
clear
close all

M = 8;
signalSize = 4000;
iter = 100;
beta = 0.9;

W = [0.4 1 -0.3];
H = [0.1 -0.2 0.5 1 -1 -0.5 0.3];

Ru = zeros(M,M);
for i=1:M
    for j=1:M
        if i==j
            Ru(i,j) = 10.56;
        elseif abs(i-j) == 1
            Ru(i,j) = -1.08;
        elseif abs(i-j) == 2
            Ru(i,j) = -5.8;
        elseif abs(i-j) == 3
            Ru(i,j) = 1.4;
        elseif abs(i-j) == 4
            Ru(i,j) = 0.6;
        elseif abs(i-j) == 5
            Ru(i,j) = -0.44;
        elseif abs(i-j) == 6
            Ru(i,j) = 0.12;
        elseif abs(i-j) == 7
            Ru(i,j) = 0;
        end
    end
end

missAdjTheory = 0.05;
miu = 2*missAdjTheory/trace(Ru);
miuHat = 2*missAdjTheory/M;

P = zeros(M, 1);
for i=1:M
    if i == 1
        P(i,1) = 0.4 * Ru(1,i) + Ru(1,i+1) - 0.3 * Ru(1,i+2);
    elseif i == 2
        P(i,1) = 0.4 * Ru(1,i) + Ru(1,i-1) - 0.3 * Ru(1,i);
    elseif i > 2
        P(i,1) = 0.4 * Ru(1,i) + Ru(i-1,1) - 0.3 * Ru(1,i-2);
    end
end

%DCT Transform:
TDCT = zeros(M,M);
for i=1:M
   if i == 1
       TDCT(1,:) = 1/sqrt(M);
   else
       for j=1:M
          TDCT(i,j) = sqrt(2/M)*cos((pi*(2*j-1)*(i-1))/(2*M)); 
       end
   end
end

%DST Transform:
TDST = zeros(M,M);
for i=1:M
   if i == 1
       TDST(1,:) = 1/sqrt(M);
   else
       for j=1:M
          TDST(i,j) = sqrt(2/(M+1))*sin((pi*(j)*(i))/(M+1)); 
       end
   end
end

%TDF Transform:
TDFT = zeros(M,M);
for i=1:M
   if i == 1
       TDFT(1,:) = 1/sqrt(M);
   else
       for j=1:M
          TDFT(i,j) = sqrt(1/M)*exp((-1i*2*pi*(j-1)*(i-1))/(M)); 
       end
   end
end

for i = 1:iter
    
   %Creating u[n]
   mSignal = 0;
   sdSignal = sqrt(4);
   input = normrnd(mSignal, sdSignal, 1, signalSize);
   u = filter(H, 1, input);
   
   %Creating uNew : output of w filter
   uNew = filter(W, 1, u);
   
   %Creating noise:
    mNoise = 0;
    sdNoise = sqrt(0.05);
    noiseSize = signalSize;
    noise = normrnd(mNoise, sdNoise, 1, noiseSize); 
    
    %Creating x: result of adding white noise to uNew
    x = uNew + noise;
    
    %Performing LMS:
    u = u.';
    x = x.';
    wLMS = zeros(M,1);
    for n=M:signalSize
       uvec = u(n:-1:n-M+1);
       e = x(n) - wLMS' * uvec;
       JLMS(n,i) = e*conj(e);
       wLMS = wLMS + miu * uvec * conj(e);
    end
    
    %Performing DCT
    wDCT = zeros(M,1);
    sigmaDCT = ones(M,M);
    for n=M:signalSize
       uvec = u(n:-1:n-M+1);
       uVec = TDCT * uvec;
       e = x(n) - wDCT' * uVec;
       for k=1:M
           wDCT(k) = wDCT(k) + miuHat * uVec(k) * conj(e) / sigmaDCT(k,k);
           sigmaDCT(k,k) = beta * sigmaDCT(k,k) + (1-beta) * (abs(uVec(k)).^2);
       end
       JDCT(n,i) = e*conj(e);
    end
    
    %Performing DST
    wDST = zeros(M,1);
    sigmaDST = ones(M,M);
    for n=M:signalSize
       uvec = u(n:-1:n-M+1);
       uVec = TDST * uvec;
       e = x(n) - wDST' * uVec;
       for k=1:M
           wDST(k) = wDST(k) + miuHat * uVec(k) * conj(e) / sigmaDST(k,k);
           sigmaDST(k,k) = beta * sigmaDST(k,k) + (1-beta) * (abs(uVec(k)).^2);
       end
       JDST(n,i) = e*conj(e);
    end
    
    %Performing DFT
    wDFT = zeros(M,1);
    sigmaDFT = ones(M,M);
    for n=M:signalSize
       uvec = u(n:-1:n-M+1);
       uVec = TDFT * uvec;
       e = x(n) - wDFT' * uVec;
       for k=1:M
           wDFT(k) = wDFT(k) + miuHat * uVec(k) * conj(e) / sigmaDFT(k,k);
           sigmaDFT(k,k) = beta * sigmaDFT(k,k) + (1-beta) * (abs(uVec(k)).^2);
       end
       JDFT(n,i) = e*conj(e);
    end
end

expJLMS = sum(JLMS,2)/iter;
expJLMS(1:M-1)=[];

expJDCT = sum(JDCT,2)/iter;
expJDCT(1:M-1)=[];

expJDST = sum(JDST,2)/iter;
expJDST(1:M-1)=[];

expJDFT = sum(JDFT,2)/iter;
expJDFT(1:M-1)=[];

Jmin = 14.426 - P' * inv(Ru) * P;
JminVec = ones(1,size(expJLMS,1))*Jmin;

temp = 0;
for i=1:100
   temp = temp + expJLMS(size(expJLMS,1)-i);  
end
JexcessLMS = temp/100 - Jmin;
misAdjLMS = (JexcessLMS/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment of LMS is :%f \n ',signalSize, iter, misAdjLMS); 

temp = 0;
for i=1:100
   temp = temp + expJDCT(size(expJDCT,1)-i);  
end
JexcessDCT = temp/100 - Jmin;
misAdjDCT = (JexcessDCT/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment of DCT is :%f \n ',signalSize, iter, misAdjDCT); 

temp = 0;
for i=1:100
   temp = temp + expJDST(size(expJDST,1)-i);  
end
JexcessDST = temp/100 - Jmin;
misAdjDST = (JexcessDST/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment of DST is :%f \n ',signalSize, iter, misAdjDST); 

temp = 0;
for i=1:100
   temp = temp + expJDFT(size(expJDFT,1)-i);  
end
JexcessDFT = temp/100 - Jmin;
misAdjDFT = (JexcessDFT/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment of DFT is :%f \n ',signalSize, iter, misAdjDFT); 

%Plotting the learning curve
figure
nExpJLMS = 0:size(expJLMS,1)-1;
plot(nExpJLMS, expJLMS, 'Color', [147/255, 112/255, 219/255]);
hold on 

nExpJDCT = 0:size(expJDCT,1)-1;
plot(nExpJDCT, expJDCT, 'Color', [102/255, 205/255, 170/255]);

nExpJDST = 0:size(expJDST,1)-1;
plot(nExpJDST, expJDST, 'Color', [106/255, 90/255, 205/255]);

nExpJDFT = 0:size(expJDFT,1)-1;
plot(nExpJDFT, expJDFT, 'Color', [255/255, 215/255, 0/255]);

nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);

title('Learning curve with optimum mu');
xlabel('n');
ylabel('E{J(w(n))}');
legend('LMS','DCT','DST','DFT','Jmin');
hold off

% finding spread value:
RuVec = TDCT * Ru * TDCT;
eignValues = eig(sigmaDCT^(-0.5)*RuVec*sigmaDCT^(-0.5));
lambdaMax = max(eignValues);
lambdaMin = min(eignValues);
spreadValue = lambdaMin / lambdaMax;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, spread value of DCT is :%f \n ',signalSize, iter, spreadValue); 

RuVec = TDST * Ru * TDST;
eignValues = eig(sigmaDST^(-0.5)*RuVec*sigmaDST^(-0.5));
lambdaMax = max(eignValues);
lambdaMin = min(eignValues);
spreadValue = lambdaMin / lambdaMax;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, spread value of DST is :%f \n ',signalSize, iter, spreadValue); 

RuVec = TDFT * Ru * TDFT;
eignValues = eig(sigmaDFT^(-0.5)*RuVec*sigmaDFT^(-0.5));
lambdaMax = max(eignValues);
lambdaMin = min(eignValues);
spreadValue = lambdaMin / lambdaMax;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, spread value of DFT is :%f \n ',signalSize, iter, spreadValue); 

