function [P,s0] = newConstrainedReduction(A,y,l,u)    
    [n m] = size(A);
    P = 1:m;
    [Q R] = qr(A);
    %flops(2/3*n^3)
    
    y = Q'*y;
    %addflops(flops_mul(Q',y));
    
    ybar=y;
    L = inv(R)';
    %addflops(n^3/3); %flops for inverse of upper-triangular
    
    s0 = zeros(1,m);
    columnMap = 1:n;
    
    for level = n:-1:2
        sig = -1;
        a = -1;
        best = -1;
        for i = 1:level
            temp = y(i:level)'*L(i:level,i);
            %addflops(flops_mul(y(i:level)',L(i:level,i)));
            
            ap = max(min(round(temp),u(i)),l(i));
            
            if(ap == l(i))
                Bp = ap + 1;
            else
                if(ap == u(i))
                    Bp = ap-1;
                else
                    Bp = ap + sign(temp-ap);
                end
            end
            %Use the lower triangular structure of L to reduce computation
            colNorm = norm(L(i:level,i));
            %addflops(flops_mul(L(i:level,i)',L(i:level,i)));
            sigp = 1./colNorm*abs(temp-Bp);
            
            if(sigp > sig)
                a = ap;
                sig = sigp;
                best = i;
            end
        end       
        s0(level) = a;
        %Swap columns of L so we can apply permutations to it in the next
        %step
        [R(1:level,1:level) y(1:level) L(1:level,1:level)] = permuTest(R(1:level,1:level),y(1:level),best,L(1:level,1:level));
        
        tempL = l(best);
        tempU = u(best);
        tempC = columnMap(best);
        for i=best:level-1
            l(i) = l(i+1);
            u(i) = u(i+1);
            columnMap(i) = columnMap(i+1);
        end
        l(level) = tempL;
        u(level) = tempU;
        columnMap(level) = tempC;
        
        P(level) = columnMap(level);
        
        y(1:level-1) = y(1:level-1) - R(1:level-1,level)*a;
        %addflops(2*(level-1));
        
    end
        P(1) = columnMap(1);
        temp = y(1)*L(1,1);
        s0(1) = max(min(round(temp),u(1)),l(1));
end