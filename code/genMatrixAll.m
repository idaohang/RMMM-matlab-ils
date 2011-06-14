function A = genMatrixAll(n,matrixCase)

%This file generates test cases
switch (matrixCase)
    case (1)
        L = tril(randn(n,n));
        D = abs(diag(rand(n,1)));
        Q = L'*D*L;
        A = D^(1/2)*L';
    case (2)
        L = tril(randn(n,n));
        D = diag((n:-1:1).^(-1));
        Q = L'*D*L;
        A = D^(1/2)*L';
    case (3)
        L = tril(randn(n,n));
        D = diag((1:n).^(-1));
        Q = L'*D*L;
        A = D^(1/2)*L';
    case (4)
        L = tril(randn(n,n));
        D = 0.1*ones(1,n);
        D(1) = 200;
        D(2) = 200;
        D(3) = 200;
        D = diag(D);
        Q = L'*D*L;
        A = D^(1/2)*L';
    case (5)
        A = randn(n,n);
        [U R] = qr(A);
        D = abs(diag(rand(n,1)));
        Q = U*D*U';
        A = U'*D^(1/2);
    case (6)
        A = randn(n,n);
        [U R] = qr(A);
        D = abs(diag(rand(n,1)));
        D(1,1) = 2^(-n/4);
        D(n,n) = 2^(n/4);
        Q = U*D*U';
        A = U'*D^(1/2);
    case (7)
        A = randn(n,n);
        Q = A'*A;
    case (8)
        A = randn(n,n);
        [U,S,V] = svd(A);
        S(1,1) = S(1,1)*10^(2);
        S(n,n) = S(n,n)*10^(-2);
        A = U*S*V';
    case (9)
        A = eye(n);
        for i=1:n
            for j =i+1:n
                A(i,j) = -0.4;%-rand/10;

            end
        end
    case (10)
        Tmp = randn(n,n);
        [U R] = qr(Tmp);
        Tmp = randn(n,n);
        [V R] = qr(Tmp);
        D = eye(n);
        for i = 1:n
            D(i,i) = 10^-(((i-1)*4)/(n-1));
        end
        A = U*D*V';
    case(11)
        A = randn(n,n);
end