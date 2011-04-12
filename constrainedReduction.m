function [P,zhat] = constrainedReduction(A,y,l,u)    
%
% [R ybar zhat lbar ubar P] = constrainedReduction(A y l u) applies the
% constrained reduction from "Solving Box-Constrained Integer Least Squares
% Problems" to the matrix A and received vector y
%
% Input arguments:
%    A ---- n by n real nonsingular matrix
%    y ---- n-dimensional real vector
%    l ---- n-dimensional lower constraint vector
%    u ---- n-dimensional upper constraint vector

% Output arguments:
%    R ---- Reduced upper triangular matrix
%    ybar ---- Reducted signal vector
%    zhat ---- constrained babai point
%    lbar ---- permuted lower constraint vector
%    ubar ---- permuted upper constraint vector
%    P ---- Permutation matrix that was applied to get R from A... this
%           must be applied to the ILS solution after it is found.
    [m n] = size(A);    
    [Q R] = qr(A);
    %flops(2*n^3)
    
    ybar = Q'*y;
    %addflops(flops_mul(Q',y));
    
    P = 1:n;
    columnMap = 1:m;
    
    zhat = zeros(n,1);
    yhat = ybar;
    lbar = l;
    ubar = u;
    %Givens = eye(n);
    
    for k = n:-1:2
        if(k < n)
            yhat(1:k) = yhat(1:k) - R(1:k,k+1)*zhat(k+1);
            %addflops(k + flops_mul(R(1:k,k+1),zhat(k+1)));
        end
        
        ck = yhat(k)/R(k,k);
        zk = max(min(round(ck),ubar(k)),lbar(k));
        if(zk == ubar(k))
            zk = zk-1;
        else
            if(zk == lbar(k))
                zk = zk+1;
            else
                zk = zk + sign(ck-zk);
            end
        end
        
        a = abs(R(k,k)*(zk-ck));
        p = k;
        Rtmp = R(1:k,1:k);
        ytmp = yhat(1:k);
        
        for j = 1:k-1
            Rp = Rtmp;
            yp = ytmp;
            
            [Rp,yp] = permu(Rp,yp,j); %swap column k and j - return to upper tri
            
            ckp = yp(k)/Rp(k,k);
            zkp = max(min(round(ckp),ubar(j)),lbar(j));
            if(zkp == ubar(j))
                zkp = zkp-1;
            else
                if(zkp == lbar(j))
                    zkp = zkp+1;
                else
                    zkp = zkp + sign(ckp-zkp);
                end
            end

            ap = abs(Rp(k,k)*(zkp-ckp));           
            
            if(ap > a)
                a = ap;
                p = j;
                R(1:k,1:k) = Rp;
                yhat(1:k) = yp;
                ck = ckp;
                %Givens = Givensp;
            end
            
        end
        

        lbark = lbar(k);
        lbar(k) = lbar(p);
        lbar(p) = lbark;
        ubark = ubar(k);
        ubar(k) = ubar(p);
        ubar(p) = ubark;

        
        P(k) = columnMap(p);
        tempMap = columnMap(k);
        columnMap(k) = columnMap(p);
        columnMap(p) = tempMap;
            
        zhat(k) = max(min(round(ck),ubar(k)),lbar(k));
    end
    P(1) = columnMap(1);
    yhat(1) = yhat(1) - R(1,2)*zhat(2);
    c1 = yhat(1)/R(1,1);
    zhat(1) = max(min(round(c1),ubar(1)),lbar(1));
end