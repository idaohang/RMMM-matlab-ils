function [R Z y] = testReduction2(A,y,alpha)
%Perform LLL, then a variant of CH where alpha gives a tradeoff between LLL
%and CH
    [m n] = size(A);
    if(checkLLL(A) == 1)
        [Q R] = qr(A);
        y = Q'*y;
        Z = eye(n,n);
    else
        [R Z y] = reduction(A,y);
    end
    [~, ~, offDiagSum,~] = checkLLL(R);
    offDiagSumOrig = offDiagSum;
    k=n;
    
    yhat = y;
    zhat = zeros(n,1);
    noCH=n+1;
    while k>=1

        if k>n
            k = n;
        end
        
        if(k<noCH)
            
            %Find the column we should permute to the kth using the CH strategy
            if(k < n)
                yhat(1:k) = yhat(1:k) - R(1:k,k+1)*zhat(k+1);
            end
            ck = yhat(k)/R(k,k);
            zk = round(ck);
            zk = zk + sign(ck-zk);
            bestDist = abs(R(k,k)*(ck-zk));
            curDist = bestDist;
            p=k;
            R2=R;
            for i = 1:k-1
                Rp=R;
                yp=yhat;
                yTemp=y;
                [Rp,yp,yTemp] = permufull(Rp,yp,yTemp,i,k); %swap column k and i - return to upper tri
                ckp = yp(k)/Rp(k,k);
                zkp = round(ckp);
                zkp = zkp + sign(ckp-zkp);
                testDist = abs(Rp(k,k)*(ckp-zkp));
                
                [~, ~, offDiagSum2,~] = checkLLL(Rp);
                measure = offDiagSum2;
                
                if(testDist > bestDist && (offDiagSum == 0 || measure <= alpha || offDiagSum2 <= offDiagSum) )
                    bestDist = testDist;
                    p = i;
                    R = Rp;
                    y=yTemp;
                    Z(:,[k,p])=Z(:,[p,k]); %record the column swap in Z
                    yhat = yp;
                    ck = ckp;
                    offDiagSum = offDiagSum2;
                end
            end

            zhat(k) = round(ck);
        end
    
        k=k-1;
    end
end