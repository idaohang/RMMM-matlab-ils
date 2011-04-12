function sigmas = getInverseColumnNorms(Q,R)
    n = length(R);
    sigmas = zeros(n,n);
    Q2 = Q';
    for k = 1:n
        for j = 1:n-k+1
            sigmas(j,n-k+1) = norm(R\Q2(:,j));
        end
        [Q,R] = qrdelete(Q,R,n-k+1,'col');
        [Q,R] = qrdelete(Q,R,n-k+1,'row');
        Q2=Q';
    end
end