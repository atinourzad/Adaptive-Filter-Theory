function [rx] = rx_finder(Rx, index)
%Since Rx is a hermity matrix for any negative index it's exact positive
%couple is used. Although all rx in here are real and conjucation has no
%effect, it's being used.
if index > 0 || index == 0
    rx = Rx(index+1);
else
    rx = conj(Rx(-index+1));
end


