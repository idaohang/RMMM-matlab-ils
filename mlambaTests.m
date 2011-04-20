runs = 100;
sz = 5:5:40;
time1 = zeros(matrixCase,size(sz,2),runs);
time2 = zeros(matrixCase,size(sz,2),runs);
permuTime = zeros(matrixCase,size(sz,2),runs);
expand1 = zeros(matrixCase,size(sz,2),runs);
expand2 = zeros(matrixCase,size(sz,2),runs);

for matrixCase = 1:7
    for n=sz;
        for i = 1:runs
            switch (matrixCase)
                case (1)
                    L = tril(randn(n,n));
                    D = diag(rand(n,1));
                    Q = L'*D*L;
                case (2)
                    L = tril(randn(n,n));
                    D = diag((n:-1:1).^(-1));
                    Q = L'*D*L;
                case (3)
                    L = tril(randn(n,n));
                    D = diag((1:n).^(-1));
                    Q = L'*D*L;
                case (4)
                    L = tril            
            tic;(randn(n,n));
                    D = 0.1*ones(1,n);
                    D(1) = 200;
                    D(2) = 200;
                    D(3) = 200;
                    D = diag(D);
                case (5)
                    A = randn(n,n);
                    [U ~] = qr(A);
                    D = diag(rand(n,1));
                    Q = U*D*U';
                case (6)
                    A = randn(n,n);
                    [U ~] = qr(A);
                    D = diag(rand(n,1));
                    D(1,1) = 2^(-n/4);
                    D(n,n) = 2^(n/4);
                case (7)
                    A = randn(n,n);
                    Q = A'*A;
            end
            ahat = 100*randn(n,1);
            y = Q*ahat;
            l = -inf*(ones(1,n));
            u = inf*(ones(1,n));
            [R1 Z1 y1] = reduction(Q,y);
            
            tic;
            [zhat1, numExpanded1] = search(R1,y1,1);
            time1(matrixCase,n,i) = toc;
            expand1(matrixCase,n,i) = numExpanded1;

            tic;
            [P z] = otherConstrainedReduction(R1,y1,l,u);
            [Q2 R2] = qr(R1(:,P));
            y2 = Q'*y1;
            permuTime(matrixCase,n,i) = toc;
            
            tic;
            [zhat2,numExpanded2] = search(R2,y2,1);
            time2(matrixCase,n,i) = toc;
            expand2(matrixCase,n,i) = numExpanded2
            
        end
    end
end