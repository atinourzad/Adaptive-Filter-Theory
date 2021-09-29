%% 

clc
clear variables
close all 

signalSize = 500;
numIter = 30;
mNoise = 0; sdNoise = 1;

e = normrnd(mNoise, sdNoise, signalSize, 1);

teta1 = zeros(numIter, signalSize); teta2 = zeros(numIter, signalSize); teta3 = zeros(numIter, signalSize);

for iter=1:numIter
    
    y = zeros(signalSize,1);
    for t=1:signalSize
        y(t,1) = 2 * sin(0.25*t) + 2 * sin(0.5 * t) + 2 * sin(0.75 * t) + e(t);
    end
    n = 3;

    yLen = size(y,1);

    %Initialization:
    tetaHatInit = [0.15; 0.35; 0.5236];
    aHatInit = [0; 0; 0];
    p = 0.1 * eye(n);
    alpha0 = 0.99; alphaInf = 0.995; alpha = 0.8;
    lambda0 = 0.8; lambdaInf = 1.01; lambda = 1.2;
    
    gamma = zeros(n,1);
    gammaPrev = zeros(n,1);
    yF = zeros(signalSize, 1);
    eF = zeros(signalSize, 1);
    eBar = zeros(signalSize, 1);
    tetaHat = zeros(n, signalSize); 
    aHat = zeros(n, signalSize);
    psi = zeros(n,1);

    for t=1:yLen
        %Step1:
        %Creating gamma
        for k=1:n-1
            if (t < k+1) && (t<(2*n-k)+1) 
                gamma(k,1) = 0;
            elseif (t < k+1) && (t>(2*n-k)) 
                gamma(k,1) = - y(t-2*n+k) + (alpha ^ (2*n-k)) * eBar(t-2*n+k);
            elseif (t > k) && (t<(2*n-k)+1) 
                gamma(k,1) = - y(t-k) + (alpha ^ k) * eBar(t-k);
            elseif (t > k) && (t>(2*n-k)) 
                gamma(k,1) = - y(t-k) - y(t-2*n+k) + (alpha ^ k) * eBar(t-k) + (alpha ^ (2*n-k)) * eBar(t-2*n+k);
            end
        end
        if t > n
           gamma(n,1) = -y(t-n) + (alpha^n) * eBar(t-n);
        else
            gamma(n,1) = 0;
        end
        %Prediction error
        if t == 1
            eHat = y(t) - gamma.' * aHatInit;
        else
            if t < 2*n+1
                eHat = y(t) - (gamma.') * aHat(:,t-1);
            else
                eHat = y(t) + y(t-2*n) - (alpha ^ (2*n)) * eBar(t-2*n) - (gamma.') * aHat(:,t-1); 
            end
        end

        %Step2:
        %Jacobian Matrix
        Jm = zeros(n,n);
        if t == 1
            Jm(1,:) = (2 * sin(tetaHatInit(:,1))).';
            for i=1:n
                Jm(2,i) = 2 * Jm(1,i) * cos(tetaHatInit(i,1)) + 2 * aHatInit(1,1) * sin(tetaHatInit(i,1));
                Jm(3,i) = 2 * Jm(2,i) * cos(tetaHatInit(i,1)) - Jm(1,i) + 2 * aHatInit(2,1) * sin(tetaHatInit(i,1));
            end
        else
            Jm(1,:) = (2 * sin(tetaHat(:,t-1))).';
            for i=1:n
                Jm(2,i) = 2 * Jm(1,i) * cos(tetaHat(i,t-1)) + 2 * aHat(1,t-1) * sin(tetaHat(i,t-1));
                Jm(3,i) = 2 * Jm(2,i) * cos(tetaHat(i,t-1)) - Jm(1,i) + 2 * aHat(2,t-1) * sin(tetaHat(i,t-1));
            end
        end
        
        %Step3:
        %Derivative vector calculation:
        %eF and yF creation
        if t == 1
            yF(t) = y(t);
            eF(t) = eBar(t);
        elseif t == 2
            yF(t) = y(t) - aHat(1,t-1) * alpha * yF(t-1);
            eF(t) = eBar(t) - aHat(1,t-1) * alpha * eF(t-1);
        elseif t == 3
            yF(t) = y(t) - aHat(1,t-1) * alpha * yF(t-1) - aHat(2,t-1) * (alpha^2) * yF(t-2);
            eF(t) = eBar(t) - aHat(1,t-1) * alpha * eF(t-1) - aHat(2,t-1) * (alpha^2) * eF(t-2);
        elseif t == 4
            yF(t) = y(t) - aHat(1,t-1) * alpha * yF(t-1) - aHat(2,t-1) * (alpha^2) * yF(t-2) - aHat(3,t-1) * (alpha^3) * yF(t-3);
            eF(t) = eBar(t) - aHat(1,t-1) * alpha * eF(t-1) - aHat(2,t-1) * (alpha^2) * eF(t-2) - aHat(3,t-1) * (alpha^3) * eF(t-3);
        elseif t == 5
            yF(t) = y(t) - aHat(1,t-1) * alpha * yF(t-1) - aHat(2,t-1) * (alpha^2) * yF(t-2) - aHat(3,t-1) * (alpha^3) * yF(t-3) - aHat(2,t-1) * (alpha^4) * yF(t-4);
            eF(t) = eBar(t) - aHat(1,t-1) * alpha * eF(t-1) - aHat(2,t-1) * (alpha^2) * eF(t-2) - aHat(3,t-1) * (alpha^3) * eF(t-3) - aHat(2,t-1) * (alpha^4) * eF(t-4);
        elseif t == 6
            yF(t) = y(t) - aHat(1,t-1) * alpha * yF(t-1) - aHat(2,t-1) * (alpha^2) * yF(t-2) - aHat(3,t-1) * (alpha^3) * yF(t-3) - aHat(2,t-1) * (alpha^4) * yF(t-4) - aHat(1,t-1) * (alpha^5) * yF(t-5);
            eF(t) = eBar(t) - aHat(1,t-1) * alpha * eF(t-1) - aHat(2,t-1) * (alpha^2) * eF(t-2) - aHat(3,t-1) * (alpha^3) * eF(t-3) - aHat(2,t-1) * (alpha^4) * eF(t-4) - aHat(1,t-1) * (alpha^5) * eF(t-5);
        else
            yF(t) = y(t) - aHat(1,t-1) * alpha * yF(t-1) - aHat(2,t-1) * (alpha^2) * yF(t-2) - aHat(3,t-1) * (alpha^3) * yF(t-3) - aHat(2,t-1) * (alpha^4) * yF(t-4) - aHat(1,t-1) * (alpha^5) * yF(t-5) - (alpha^6) * yF(t-6);
            eF(t) = eBar(t) - aHat(1,t-1) * alpha * eF(t-1) - aHat(2,t-1) * (alpha^2) * eF(t-2) - aHat(3,t-1) * (alpha^3) * eF(t-3) - aHat(2,t-1) * (alpha^4) * eF(t-4) - aHat(1,t-1) * (alpha^5) * eF(t-5) - (alpha^6) * eF(t-6);
        end

        %Gamma creation:
        if t < 3
            gammaPrev(:,1) = 0;
        elseif t == 3
            gammaPrev(:,1) = 0;
            gammaPrev(1,1) = - yF(t-2) + alphaPrev * eF(t-2);
        elseif t == 4
            gammaPrev(1,1) = - yF(t-2) + alphaPrev * eF(t-2);
            gammaPrev(2,1) = - yF(t-3) + (alphaPrev^2) * eF(t-3);
            gammaPrev(3,1) = 0;
        elseif t == 5
            gammaPrev(1,1) = - yF(t-2) + alphaPrev * eF(t-2);
            gammaPrev(2,1) = - yF(t-3) + (alphaPrev^2) * eF(t-3);
            gammaPrev(3,1) = - yF(t-4) + (alphaPrev^3) * eF(t-4);
        elseif t == 6
            gammaPrev(1,1) = - yF(t-2) + alphaPrev * eF(t-2);
            gammaPrev(2,1) = - yF(t-3) + (alphaPrev^2) * eF(t-3) - yF(t-5) + (alphaPrev^4) * eF(t-5);
            gammaPrev(3,1) = - yF(t-4) + (alphaPrev^3) * eF(t-4);
        else
            gammaPrev(1,1) = - yF(t-2) + alphaPrev * eF(t-2) - yF(t-6) + (alphaPrev^5) * eF(t-6);
            gammaPrev(2,1) = - yF(t-3) + (alphaPrev^2) * eF(t-3) - yF(t-5) + (alphaPrev^4) * eF(t-5);
            gammaPrev(3,1) = - yF(t-4) + (alphaPrev^3) * eF(t-4); 
        end

        psi = Jm * gammaPrev;

        %Step4
        %Parameter updating
        k = (p * psi) / (lambda + (psi.') * p * psi);
        p = (p - k * (psi.') * p) / lambda;
        if t==1
            tetaHat(:,t) = tetaHatInit(:,1) + k * eHat;
        else
            tetaHat(:,t) = tetaHat(:,t-1) + k * eHat;
        end

        %Step5:
        %Conversion of teta to a
        aHat(1,t) = -2 * sum(cos(tetaHat(:,t)),1);
        aHat(2,t) = 3 + 4 * cos(tetaHat(3,t)) * (cos(tetaHat(1,t))+cos(tetaHat(2,t))) + 4 * cos(tetaHat(1,t)) * cos(tetaHat(2,t));
        aHat(3,t) = -4 * cos(tetaHat(1,t)) - 4 * cos(tetaHat(2,t)) - 4 * cos(tetaHat(3,t)) * (1+2*cos(tetaHat(1,t))*cos(tetaHat(2,t)));

        %Step6:
        %Posteriori prediction error
        if t < 2*n+1
            eBar(t) = y(t) - gamma.' * aHat(:,t);
        else
            eBar(t) = y(t) + y(t-2*n) - (alpha ^ (2*n)) * eBar(t-2*n) - gamma.' * aHat(:,t); 
        end

        %Step7:
        %Constant updaing
        alphaPrev = alpha;
        alpha = alphaInf - (alphaInf - alpha) * alpha0;
        lambda = lambdaInf - (lambdaInf - lambda) * lambda0;
    end

    teta1(iter,:) = abs(tetaHat(1,:));
    teta2(iter,:) = abs(tetaHat(2,:));
    teta3(iter,:) = abs(tetaHat(3,:));  
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

title('Frequency Evolition using CYL''s Algorithm');
xlabel('Iteration Number');
ylabel('Frequency');

hold off

