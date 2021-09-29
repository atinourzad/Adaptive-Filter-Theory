function [x_ar5] = AR5_finder(v, W, sd_noise, noise_size)
%For AR(5): x[n] = W(1)*x[n-1] + W(2)*x[n-2] + ... + W(5)*x[n-5] + s_d_v*v[n].
x_ar5 = zeros(noise_size, 1);
for i=1:noise_size
    if i==1 %Since a negative index has no meaning and it's hypothesized that signal has value equals to zero there, these conditions are applied. 
        x_ar5(i) = sd_noise * v(i);
    elseif i==2
        x_ar5(i) = W(1) * x_ar5(i-1) + sd_noise * v(i);
    elseif i==3
        x_ar5(i) = W(1) * x_ar5(i-1) + W(2) * x_ar5(i-2) + sd_noise * v(i);
    elseif i==4
        x_ar5(i) = W(1) * x_ar5(i-1) + W(2) * x_ar5(i-2) + W(3) * x_ar5(i-3) + sd_noise * v(i);
    elseif i==5
        x_ar5(i) = W(1) * x_ar5(i-1) + W(2) * x_ar5(i-2) + W(3) * x_ar5(i-3) + W(4) * x_ar5(i-4) + sd_noise * v(i);
    else
        x_ar5(i) = W(1) * x_ar5(i-1) + W(2) * x_ar5(i-2) + W(3) * x_ar5(i-3) + W(4) * x_ar5(i-4) + W(5) * x_ar5(i-5) + sd_noise * v(i);
    end
end
end

