function [Rx_ar] = AR_Rx_finder(Rx, size)
%Making Rx for an AR Process using rx values. Since for any AR(M),
%rx(-M+1), rx(-M),..., rx(-1),rx(0),...,rx(M-1) is needed, Rx is created
%using them.
Rx_ar = zeros(size);
for i=1:size
    for j=1:size
        Rx_ar(i,j) = rx_finder(Rx, j-i);
    end
end 
end

