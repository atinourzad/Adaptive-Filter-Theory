%% B
clc
clear
close all

%White noise generation with mean 0 and s.d 1
noise_size = 100;
v = normrnd(0, 1, noise_size, 1);

%MA process: x[n] = v[n] + 0.1v[n-1] + 0.25v[n-2]
x = zeros(noise_size, 1);
for i=1:noise_size
    if i == 1
        x(i) = v(i);
    elseif i == 2
        x(i) = v(i) + 0.1 * v(i-1);
    else
        x(i) = v(i) + 0.1 * v(i-1) + 0.25 * v(i-2);
    end
end

%AR(2) Process: x[n] = 0.0906x[n-1] + 0.2225x[n-2] + 1.00277v[n]
W = [0.0906; 0.2225];
sd_noise = 1.00277;
x_ar2 = AR2_finder(v, W, sd_noise, noise_size);

%AR(5) Process: x[n] = 0.1001x[n-1] + 0.2393x[n-2] - 0.0493x(n-3) - 0.0521x(n-4) + 0.0176x(n-5) + 1.00008v[n]
W = [0.1001;0.2393;-0.0493;-0.0521;0.0176];
sd_noise = 1.00008;
x_ar5 = AR5_finder(v, W, sd_noise, noise_size);

%AR(10): x[n] = 0.1x[n-1] + 0.24x[n-2] - 0.049x(n-3) - 0.0551x(n-4) + 0.0178x(n-5) + 0.012x(n-6) - 0.0056x[n-7] - 0.0024x[n-8] + 00015x[n-9] + 0.0004x[n-10] + v[n]
W = [0.1;0.24;-0.049;-0.0551;0.0178;0.012;-0.0056;-0.0024;00015;0.0004];
x_ar10 = AR5_finder(v, W, 1, noise_size);

figure
plot(x, 'Color', [147/255, 112/255, 219/255])
hold on 
plot(x_ar2, 'Color', [218/255, 112/255, 214/255]);
plot(x_ar5, 'Color', [102/255, 205/255, 170/255]);
plot(x_ar10, 'Color', [255/255, 215/255, 0/255]);
legend('MA Process','AR(2)','AR(5)','AR(10)');
hold off

%% C 
%Finding the mean square error
clc
close all

err_ar2 = immse(x, x_ar2); %Between x[n] and it's AR(2) approximation
err_ar5 = immse(x, x_ar5); %Between x[n] and it's AR(5) approximation
err_ar10 = immse(x, x_ar10); %Between x[n] and it's AR(10) approximation

%% D

clc
clear x_ar2 x_ar5 x_ar10
close all

%Finding Rx based on the approximation asked in question.
Rx = zeros(1,noise_size);
for i=0:noise_size-1
    Rx(1,i+1) = rx_approx(noise_size, i, x);
end
%Since Rx is a hermity matrix, i.e. for every m: rx(m) = rx*(-m), finding the
%other components are not necessary.

%AR(2) Process: 
Rx_ar2 =  AR_Rx_finder(Rx, 2); 
rx_ar2 = Rx(2:3);
W = Rx_ar2\rx_ar2';
sd_noise = sqrt(noise_sd_finder(Rx, 2, W));
x_ar2 =  AR2_finder(v, W, sd_noise, noise_size);

%AR(5) Process: 
Rx_ar5 =  AR_Rx_finder(Rx, 5); 
rx_ar5 = Rx(2:6);
W = Rx_ar5\rx_ar5';
sd_noise = sqrt(noise_sd_finder(Rx, 5, W));
x_ar5 =  AR5_finder(v, W, sd_noise, noise_size);

%AR(10) Process: 
Rx_ar10 =  AR_Rx_finder(Rx, 10); 
rx_ar10 = Rx(2:11);
W = Rx_ar10\rx_ar10';
sd_noise = sqrt(noise_sd_finder(Rx, 10, W));
x_ar10 =  AR10_finder(v, W, sd_noise, noise_size);

figure
plot(x, 'Color', [147/255, 112/255, 219/255])
hold on 
plot(x_ar2, 'Color', [218/255, 112/255, 214/255]);
plot(x_ar5, 'Color', [102/255, 205/255, 170/255]);
plot(x_ar10, 'Color', [255/255, 215/255, 0/255]);
legend('MA Process', 'AR(2)', 'AR(5)', 'AR(10)');
hold off


err_ar2 = immse(x, x_ar2); %Between x[n] and it's AR(2) approximation
err_ar5 = immse(x, x_ar5); %Between x[n] and it's AR(5) approximation
err_ar10 = immse(x, x_ar10); %Between x[n] and it's AR(10) approximation

%% D with 1000 iterations : 

clc
clear
close all

%White noise generation with mean 0 and s.d 1
noise_size = 1000;
v = normrnd(0, 1, noise_size, 1);

%MA process: x[n] = v[n] + 0.1v[n-1] + 0.25v[n-2]
x = zeros(noise_size, 1);
for i=1:noise_size
    if i == 1
        x(i) = v(i);
    elseif i == 2
        x(i) = v(i) + 0.1 * v(i-1);
    else
        x(i) = v(i) + 0.1 * v(i-1) + 0.25 * v(i-2);
    end
end

%Finding Rx based on the approximation asked in question.
Rx = zeros(1,noise_size);
for i=0:noise_size-1
    Rx(1,i+1) = rx_approx(noise_size, i, x);
end
%Since Rx is a hermity matrix, i.e. for every m: rx(m) = rx*(-m), finding the
%other components are not necessary.

%AR(2) Process: 
Rx_ar2 =  AR_Rx_finder(Rx, 2); 
rx_ar2 = Rx(2:3);
W = Rx_ar2\rx_ar2';
sd_noise = sqrt(noise_sd_finder(Rx, 2, W));
x_ar2 =  AR2_finder(v, W, sd_noise, noise_size);

%AR(5) Process: 
Rx_ar5 =  AR_Rx_finder(Rx, 5); 
rx_ar5 = Rx(2:6);
W = Rx_ar5\rx_ar5';
sd_noise = sqrt(noise_sd_finder(Rx, 5, W));
x_ar5 =  AR5_finder(v, W, sd_noise, noise_size);

%AR(10) Process: 
Rx_ar10 =  AR_Rx_finder(Rx, 10); 
rx_ar10 = Rx(2:11);
W = Rx_ar10\rx_ar10';
sd_noise = sqrt(noise_sd_finder(Rx, 10, W));
x_ar10 =  AR10_finder(v, W, sd_noise, noise_size);

figure
plot(x, 'Color', [147/255, 112/255, 219/255])
hold on 
plot(x_ar2, 'Color', [218/255, 112/255, 214/255]);
plot(x_ar5, 'Color', [102/255, 205/255, 170/255]);
plot(x_ar10, 'Color', [255/255, 215/255, 0/255]);
legend('MA Process', 'AR(2)', 'AR(5)', 'AR(10)');
hold off


err_ar2 = immse(x, x_ar2); %Between x[n] and it's AR(2) approximation
err_ar5 = immse(x, x_ar5); %Between x[n] and it's AR(5) approximation
err_ar10 = immse(x, x_ar10); %Between x[n] and it's AR(10) approximation

