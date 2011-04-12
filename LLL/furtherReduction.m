clc; clear;

n=5;
delta = 1;
A = rand(n,n)*100;

[U,S,V] = svd(A);
S(1,1) = S(1,1)*10^(2);
S(n,n) = S(n,n)*10^(-2);
A = U*S*V';

[Q R] = qr(A);
R1 = R;
Z = eye(n,n);
k=n;

while k>=2
    
    if k>=n
        k = n;
    end
    
    % Size reduction
    [R,Z] = IGT(R,Z,k);
    
    % Column swap
    gamma = R(k,k)^2 + R(k-1,k)^2;
    if gamma < delta*R(k-1,k-1)^2
        [R,Z] = swap(R,k-1,Z);
        k = k + 1;
    else    
        k = k - 1;
    end 
    
end

