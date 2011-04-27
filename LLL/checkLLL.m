%Returns 1 if R is LLL reduced, 0 otherwise
function [check,row,offdiag,avgAngle] = checkLLL(R)
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
        for j = i+1:n
            if(abs(R(i,j)) > 0.5*abs(R(i,i)))
                row = row + abs(R(i,j)) - 0.5*abs(R(i,i));
                check=0;
            end
        end
    end
    
    %Compute the average angle between columns... is invariant under
    %permutations and works out to a pretty consistent average for all
    %problem sizes surprisingly
    sum=0;
    count=0;
    for i = 1:n
        for j = i+1:n
            count=count+1;
            sum = sum + acos((R(:,i)./norm(R(:,i)))'*(R(:,j)./norm(R(:,j))));
        end
    end
    avgAngle = sum/count;
    
end