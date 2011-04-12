function [eVals,F] = preProcessLB(R)
    %First compute inverse of R
    m = length(R);
    F = inv(R);
    eVals = zeros(m,1);
    for i=1:m
       eVals(i) = min(eig(R(1:i,1:i)'*R(1:i,1:i)));
    end
end