%% 

clc
clear
close all 

signalSize = 500;
numIter = 30;
mNoise = 0; sdNoise = 1;

y = zeros(signalSize,1);

teta1 = zeros(numIter, signalSize); teta2 = zeros(numIter, signalSize); teta3 = zeros(numIter, signalSize);

for iter=1:numIter
    e = normrnd(mNoise, sdNoise, signalSize, 1);
    for t=1:signalSize
        y(t,1) = 2 * sin(0.5*t) + 2 * sin(t) + 2 * sin(2 * t) + 2 * sin(3 * t) + e(t);
    end
    n = 3;
    yLen = size(y,1);

    %Initialization:
    tetaHatInit = [2.6180; 1.5708; 0.5236];
    p = 0.1 * eye(n);
    alpha0 = 0.99; alphaInf = 0.995; alpha = 0.8;
    lambda0 = 0.9; lambdaInf = 1; lambda = 0.98;
    beta = 0.9999;

    z = zeros(2,n);
    wAlpha = zeros(2,n); wBeta = zeros(2,n);
    tetaHat = zeros(n, signalSize); 
    psi = zeros(n,1);

    %Main Loop:
    for t=1:signalSize

        %Step1:
        %Prediction error
        x = y(t);
        for i=1:n
            if t == 1
                x = [2*(alpha - beta)* cos(tetaHatInit(i,1)) beta^2-alpha^2] * z(:,i) + x;
            else
                x = [2*(alpha - beta)* cos(tetaHat(i,t-1)) beta^2-alpha^2] * z(:,i) + x;
            end
        end
        eHat = x;

        %Step2:
        %Computing derivative:
        for i=1:n
           eFiBeta = [beta 0] * wBeta(:, i);
           eFiAlpha = [alpha 0] * wAlpha(:, i);
           if t==1
               psi(i,1) = -2 * (eFiBeta - eFiAlpha) * sin(tetaHatInit(i,1));
           else
               psi(i,1) = -2 * (eFiBeta - eFiAlpha) * sin(tetaHat(i,t-1));
           end
        end
        %Step3:
        %Teta updating
        k = (p * psi) / (lambda + (psi.') * p * psi);
        p = (p - k * (psi.') * p) / lambda;
        if t==1
            tetaHat(:,t) = tetaHatInit(:,1) + k * eHat;
        else
            tetaHat(:,t) = tetaHat(:,t-1) + k * eHat;
        end

        %Step4:
        %Posteriori prediction error
        x = y(t);
        for i=1:n
            dummyX = x;
            x = [2*(alpha - beta)* cos(tetaHat(i,t)),beta^2-alpha^2] * z(:,i) + x;
            z(:,i) = [2*alpha*cos(tetaHat(i,t)),-alpha^2;1,0] * z(:,i) + [1;0] * dummyX;
        end
        eBar = x;
        for i=1:n
           wAlpha(:,i) = [2*alpha*cos(tetaHat(i,t)),-alpha^2;1,0] * wAlpha(:,i) + [1;0] * eBar;
           wBeta(:,i) = [2*alpha*cos(tetaHat(i,t)),-alpha^2;1,0] * wBeta(:,i) + [1;0] * eBar;
        end

        alpha = alphaInf - (alphaInf - alpha) * alpha0;
        lambda = lambdaInf - (lambdaInf - lambda) * lambda0;
    end
    
    for i=1:n
       if abs(tetaHat(i,signalSize)) < 0.8
            teta1(iter,:) = abs(tetaHat(i,:));
       elseif abs(tetaHat(i,signalSize)) < 1.6
           teta2(iter,:) = abs(tetaHat(i,:));
       else
           teta3(iter,:) = abs(tetaHat(i,:));
       end
    end    
end

expTeta1 = sum(teta1,1)/numIter;
expTeta2 = sum(teta2,1)/numIter;
expTeta3 = sum(teta3,1)/(numIter);


figure

nExpTeta1 = 0:size(expTeta1,2)-1;
plot(nExpTeta1, expTeta1, 'Color', [147/255, 112/255, 219/255]);

hold on 
nExpTeta2 = 0:size(expTeta2,2)-1;
plot(nExpTeta2, expTeta2, 'Color', [102/255, 205/255, 170/255]);

nExpTeta3 = 0:size(expTeta3,2)-1;
plot(nExpTeta3, expTeta3, 'Color', [255/255, 0/255, 255/255]);

title('Frequency Evolition using GL Algorithm');
xlabel('Iteration Number');
ylabel('Frequency');

hold off