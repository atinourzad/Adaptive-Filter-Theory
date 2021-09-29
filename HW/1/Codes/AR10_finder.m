function [x_ar10] = AR10_finder(v, W, sd_noise, noise_size)
%For AR(10): x[n] = W(1)*x[n-1] + W(2)*x[n-2] + ... + W(10)*x[n-10] + s_d_v*v[n].
x_ar10 = zeros(noise_size, 1);
for i=1:noise_size
    if i==1 %Since a negative index has no meaning and it's hypothesized that signal has value equals to zero there, these conditions are applied. 
        x_ar10(i) = sd_noise * v(i);
    elseif i==2
        x_ar10(i) = W(1) * x_ar10(i-1) + sd_noise * v(i);
    elseif i==3
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + sd_noise * v(i);
    elseif i==4
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + sd_noise * v(i);
    elseif i==5
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + W(4) * x_ar10(i-4) + sd_noise * v(i);
    elseif i==6
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + W(4) * x_ar10(i-4) + W(5) * x_ar10(i-5) + sd_noise * v(i);
    elseif i==7
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + W(4) * x_ar10(i-4) + W(5) * x_ar10(i-5) + W(6) * x_ar10(i-6) + sd_noise * v(i);
    elseif i==8
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + W(4) * x_ar10(i-4) + W(5) * x_ar10(i-5) + W(6) * x_ar10(i-6) + W(7) * x_ar10(i-7) + sd_noise * v(i);
    elseif i==9
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + W(4) * x_ar10(i-4) + W(5) * x_ar10(i-5) + W(6) * x_ar10(i-6) + W(7) * x_ar10(i-7) + W(8) * x_ar10(i-8)  + sd_noise * v(i);
    elseif i==10
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + W(4) * x_ar10(i-4) + W(5) * x_ar10(i-5) + W(6) * x_ar10(i-6) + W(7) * x_ar10(i-7) + W(8) * x_ar10(i-8) + W(9) * x_ar10(i-9) + sd_noise * v(i);
    else
        x_ar10(i) = W(1) * x_ar10(i-1) + W(2) * x_ar10(i-2) + W(3) * x_ar10(i-3) + W(4) * x_ar10(i-4) + W(5) * x_ar10(i-5) + W(6) * x_ar10(i-6) + W(7) * x_ar10(i-7) + W(8) * x_ar10(i-8) + W(9) * x_ar10(i-9) + W(10) * x_ar10(i-10) + sd_noise * v(i);
    end
end

end

