function [R,Z,y] = reduction(B,y)
%
% [R,Z,y] = reduction(B,y) computes the LLL-QRZ factorization:
%           Q'*B*Z = [R; 0] and computes Q'*y. The orthogonal matrix Q 
%           is not produced. Its goal is to reduce a general integer
%           least squares problem to an upper triangular one.
%
% Input arguments:
%    B --- m by n real matrix with full column rank
%    y --- m-dimensional real vector to be transformed to Q'y
%
% Output arguments:
%    R --- n by n LLL-reduced upper triangular matrix
%    Z --- n by n unimodular matrix, i.e., an integer matrix with |det(Z)|=1
%    y --- m-vector transformed from the input y by Q', i.e., y := Q'*y

% Copyright (c) 2006. Xiao-Wen Chang and Tianyang Zhou
% Version 1.0, October 2006.


% Check input arguments
if nargin < 2  % input error
    error('Not enough arguments!')
end

[m,n] = size(B);
if m < n  % input error
    error('Input matrix is column-rank deficient!')
end

[m2,n2] = size(y);
if m ~= m2 | n2 ~= 1  % input error
    error('Input arguments have a dimension error!')
end

% QR with minimum-column pivoting
[R,y,piv] = qrmcp(B,y);

% Initialize the unimodual matrix Z
Z = eye(n);

% Update Z using piv
for j = 1:n
    Z(:,[j,piv(j)]) = Z(:,[piv(j),j]);
end

while 1
    
    % Find R(:,j1:j2) to do off-diagonal size reduction and column reordering
    
    % Find j2
    j2 = 0;
    minratio = 1;
    for j = 2:n
        if abs(R(j,j)/R(j-1,j-1)) < minratio
            i = j - 1;
            alpha = round(R(i,j)/R(i,i));
            eta = R(i,j) - alpha*R(i,i);
            ratio = sqrt(R(j,j)^2 + eta^2)/abs(R(i,i));
            if ratio < minratio
                j2 = j;
                minratio = ratio;
            end
        end  
    end
          
    % If j2 is not found, exit
    if j2 == 0
        break;
    end
                
    % Perform off-diagonal size reduction for R(:,j2)
    for i = j2-1:-1:1
        alpha = round(R(i,j2)/R(i,i));
        if alpha ~= 0
            R(1:i,j2) = R(1:i,j2) - alpha*R(1:i,i);
            Z(:,j2) = Z(:,j2) - alpha*Z(:,i);
        end
    end

    % Find j1
    j1 = j2 - 1;
    if j1 > 1
        rho = norm(R(j1-1:j2,j2))^2;
        while rho < R(j1-1,j1-1)^2
            j1 = j1-1;
            if j1 == 1
                break;
            end
            rho = rho + R(j1-1,j2)^2;
        end
    end

    % Perform off-diagonal size reduction for R(:,j1+1:j2-1)
    for j = j1+1:j2-1
        for i = j-1:-1:1
            alpha = round(R(i,j)/R(i,i));
            if alpha ~= 0
                R(1:i,j) = R(1:i,j) - alpha*R(1:i,i);
                Z(:,j) = Z(:,j) - alpha*Z(:,i);
            end
        end
    end

    % Perform QR factorization and do column reordering for R(:,j1:j2)
    [R(j1:j2,j1:j2),T,piv] = qrmcp(R(j1:j2,j1:j2),[y(j1:j2),R(j1:j2,j2+1:n)]);
    y(j1:j2) = T(:,1);
    R(j1:j2,j2+1:n) = T(:,2:n-j2+1);
    
    % Reorder columns of R(1:j1-1,j1:j2) and Z using piv
    for j = j1:j2
        k = j1 + piv(j-j1+1) - 1;
        R(1:j1-1,[j,k]) = R(1:j1-1,[k,j]);
        Z(:,[j,k]) = Z(:,[k,j]);
    end

    % Perform off-diagonal size reduction for R(:,j1:j2)
    for j = j1:j2
        for i = j-1:-1:1
            alpha = round(R(i,j)/R(i,i));
            if alpha ~= 0
                R(1:i,j) = R(1:i,j) - alpha*R(1:i,i);
                Z(:,j) = Z(:,j) - alpha*Z(:,i);
            end
        end
    end 
        
end

% Perform off-diagonal size reduction for R
for j = 2:n
    for i = j-1:-1:1
        alpha = round(R(i,j)/R(i,i));
        if alpha ~= 0
            R(1:i,j) = R(1:i,j) - alpha*R(1:i,i);
            Z(:,j) = Z(:,j) - alpha*Z(:,i);
        end
    end
end
