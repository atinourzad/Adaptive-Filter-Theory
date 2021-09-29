function [x_ar2] = AR2_finder(v, W, sd_noise, noise_size)
%For AR(2): x[n] = W(1)*x[n-1] + W(2)*x[n-2] + s_d_v*v[n].
x_ar2 = zeros(noise_size, 1);
for i=1:noise_size
    if i==1 %Since a negative index has no meaning and it's hypothesized that signal has value equals to zero there, these conditions are applied. 
        x_ar2(i) = sd_noise * v(i);
    elseif i==2
        x_ar2(i) = W(1) * x_ar2(i-1) + sd_noise * v(i);
    else
        x_ar2(i) = W(1) * x_ar2(i-1) + W(2) * x_ar2(i-2) + sd_noise * v(i);
    end
end
end

