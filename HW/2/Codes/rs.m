function [rs_n] = rs(n)
rs_n = 25*( (1/3)*power(0.5, abs(n)) - (5/21)*power(0.4, abs(n)) );
end

