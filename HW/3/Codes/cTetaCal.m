function [cTeta] = cTetaCal(teta, M)
cTeta = zeros(M, size(teta,1));
for i=1:size(teta,1)
    for j=0:M-1
        cTeta(j+1,i) = exp(1i*pi*cos(teta(i)*pi/180)*j); 
    end
end
end

