%% Part 3

clc
clear 
close all

M = 35;
signalSize = 2000;
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

Jmin = 1 - P' * inv(Rx) * P;

% Miss adjustment = %1 
expJRLS1 = myRLS(M, signalSize, iter, 1/100);

temp = 0;
for i=1:100
   temp = temp + expJRLS1(size(expJRLS1,1)-i);  
end

Jexcess1 = temp/100 - Jmin;
misadj1 = (Jexcess1/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f \n ',signalSize, iter, misadj1); 

% Miss adjustment = %5
expJRLS2 = myRLS(M, signalSize, iter, 5/100);

temp = 0;
for i=1:100
   temp = temp + expJRLS2(size(expJRLS2,1)-i);  
end

Jexcess2 = temp/100 - Jmin;
misadj2 = (Jexcess2/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f \n ',signalSize, iter, misadj2); 

% Miss adjustment = %20
expJRLS3 = myRLS(M, signalSize, iter, 20/100);

temp = 0;
for i=1:100
   temp = temp + expJRLS3(size(expJRLS3,1)-i);  
end

Jexcess3 = temp/100 - Jmin;
misadj3 = (Jexcess3/Jmin)*100;
fprintf('With %d number of iterations (signal size) and %d number of re-dos, misadjustment is :%f \n ',signalSize, iter, misadj3); 

JminVec = ones(1,size(expJRLS1,1))*Jmin;

%Plotting the learning curve
figure
nExpJ1 = 0:size(expJRLS1,1)-1;
plot(nExpJ1, expJRLS1, 'Color', [147/255, 112/255, 219/255]);
hold on 
nExpJ2 = 0:size(expJRLS2,1)-1;
plot(nExpJ2, expJRLS2, 'Color', [102/255, 205/255, 170/255]);

nExpJ3 = 0:size(expJRLS3,1)-1;
plot(nExpJ3, expJRLS3, 'Color', [255/255, 0/255, 255/255]);

nJmin = 0:size(JminVec,2)-1;
plot(nJmin, JminVec, 'Color', [255/255, 215/255, 0/255]);

title('Learning curve');
xlabel('n');
ylabel('E{J(w(n))}');
legend('RLS 1%','RLS 5%','RLS 20%','Jmin');
hold off


