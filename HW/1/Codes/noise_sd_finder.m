function [sd_noise] = noise_sd_finder(Rx, N, W)
%Based on standard deviation of v[n] is defined using W and rx, it's calculated. 
sd_noise = Rx(1);
for i = 1:N
    sd_noise = sd_noise - Rx(i+1) * W(i);
end

