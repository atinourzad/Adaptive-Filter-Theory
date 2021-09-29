function [rx] = rx_approx(N, m, x)
%Using correlation function estination.
rx = 0;
for k = 1:N-m
    rx = rx + (1/(N-m)) * ( x(m+k) * conj(x(k)) );
end
