function [R Z y] = newReduction(A,y)
    delta=1;
    [m n] = size(A);
    [Q R] = qr(A);
    y = Q'*y;
    Z = eye(n,n);
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
            p=k;
            for i = 1:k-1
                Rp=R;
                yp=yhat;
                yTemp=y;
                [Rp,yp,yTemp] = permufull(Rp,yp,yTemp,i,k); %swap column k and i - return to upper tri
                ckp = yp(k)/Rp(k,k);
                zkp = round(ckp);
                zkp = zkp + sign(ckp-zkp);
                testDist = abs(Rp(k,k)*(ckp-zkp));
                if(testDist > bestDist)
                    bestDist = testDist;
                    p = i;
                    R = Rp;
                    y=yTemp;
                    Z(:,[k,p])=Z(:,[p,k]); %record the column swap in Z
                    yhat = yp;
                    ck = ckp;
                end
            end
            zhat(k) = round(ck);
        end
        
        if(k < n)
            %Size reduce previous column
            [R Z] = IGT(R,Z,k+1);
            
            % Column swap previous if necessary
            gamma = R(k+1,k+1)^2 + R(k,k+1)^2;
            if gamma < delta*R(k,k)^2
                [R,Z,y] = swap(R,k,Z,y);
                zhat([k+1,k]) = zhat([k,k+1]);
                %After swap, size reduce again
                [R Z] = IGT(R,Z,k+1);
                noCH = min(noCH,k);
                k=k+1;
            else
                k=k-1;
            end
        else
            k=k-1;
        end        
    end
    for i = 1:n
        [R Z] = IGT(R,Z,i);
    end
end