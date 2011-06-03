function [P,s0] = otherConstrainedReduction(A,y,l,u)
    %flops(0);
    [n m] = size(A);
    I=eye(n);
    P = 1:m;
    Index = 1:m;
    G = pinv(A)';
    %addflops(2*n^3);

    s0 = zeros(1,m);
    for L = m:-1:1

        sig = -1;
        a = -1;
        best = -1;
        bestColNorm = 0;
        count=0;
        for i = Index
            count=count+1;
            temp = y'*G(:,i);
            %addflops(flops_mul(y',G(:,i)));

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
            colNorm = norm(G(:,i));
            %addflops(flops_mul(G(:,i)',G(:,i)));
            sigp = 1/colNorm*abs(temp-Bp);
            fprintf('%i, %i, %f, %f, %i\n',L,i,sigp,temp,Bp);
            
            if(sigp > sig)
                a = ap;
                sig = sigp;
                best = i;
                bestPos = count;
                bestColNorm = colNorm*colNorm;
            end
        end
        
        P(L) = best;
        s0(L) = a;
        
        Index(bestPos) = [];       
        shiftY = (y-A(:,best)*a);
        %addflops(2*n);
        divCol = G(:,best)./bestColNorm;
        %addflops(flops_div*n);
        y = shiftY; %- divCol*(G(:,best)'*shiftY);
        %addflops(2*n + flops_mul(G(:,best)',shiftY))
        
        for i = Index
            i
            G(:,i) = G(:,i) - divCol*(G(:,best)'*G(:,i));
            %addflops(2*n + flops_mul(G(:,best)',G(:,i)));
        end
        G
        y
    end
end