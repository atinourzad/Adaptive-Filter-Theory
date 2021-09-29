function [JNew] = JCalc(varD, P, W, R, JOld)
JNew = JOld;
JNew(1, size(JOld, 2)+1) = varD - P' * W - W' * P + W' * R * W;
end

