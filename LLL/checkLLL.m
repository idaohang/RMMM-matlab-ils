%Returns 1 if R is LLL reduced, 0 otherwise
function [check,row,offdiag] = checkLLL(R)
    check = 1;
    offdiag=0;
    row=0;
    if(triu(R) ~= R)
        error('Matrix must be upper triangular');
    end
    [m n] = size(R);
    for i = 1:m
        if ( i > 1 && R(i,i)^2 + R(i-1,i)^2 < R(i-1,i-1)^2)
            offdiag = offdiag + R(i-1,i-1)^2 - (R(i,i)^2 + R(i-1,i)^2);
            check=0;
        end
        for j = 2:n
            if(i~=j && abs(R(i,j)) > 0.5*abs(R(i,i)))
                row = row + abs(R(i,j)) - 0.5*abs(R(i,i));
                check=0;
            end
        end
    end
end