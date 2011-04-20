function [R,F,piv] = qrmcp(C,F)
%
% [R,F,piv] = qrmcp(C,F) computes the QR factorization with 
%             minimum-column pivoting: Q'CP = [R; 0] and computes Q'*F. 
%             The orthogonal matrix Q is not produced. 
%
% Input arguments:
%    C --- m by n real matrix to be factorized
%    F --- m by p real matrix to be transformed to Q'F
%
% Output arguments:
%    R --- n by n real upper triangular matrix
%    F --- m by p matrix transormed from the input F by Q', i.e., F := Q'*F
%    piv - n-vector storing the information of the permutation matrix P

% Copyright (c) 2006. Xiao-Wen Chang and Tianyang Zhou
% Version 1.0,  October 2006.


% Check input arguments
[m,n] = size(C);
if m < n  % input error
    error('Matrix to be factorized is column-rank deficient!')
end;

if nargin < 2
    p = 0;
    F = zeros(m,0);
else
    [m2,p] = size(F);
    if m ~= m2  % input error
        error('Dimensions do not match!')
    end;
end

% Initialization
s = zeros(n,1);
piv = zeros(n,1);

for j = 1:n
    s(j) = norm(C(:,j))^2; % C(:,j)'*C(:,j)
end

for k = 1:n

    % Find the column with minimum 2-norm 
    [minnorm, i] = min(s(k:n));
    q = i+k-1;
    piv(k) = q;
    
    % Column interchange
    if q > k
        C(:,[k,q]) = C(:,[q,k]);
        s([k,q]) = s([q,k]);
    end

    % Compute and apply the Householder transformation  I-tau*v*v'
    if C(k,k)^2 < s(k),  % otherwise no Householder transformation is needed
        l = m+1-k;
        v = C(k:m,k);
        rho = sqrt(s(k));    % R(k,k)=norm(C(k:n,k))
        if (v(1) >= 0) 
            rho = -rho;
        end
        v(1) = v(1) - rho;   % C(k,k)+sgn(C(k,k))*norm(C(k:n,k))
        tao = -1/(rho*v(1));
        C(k,k) = rho;
        C(k:m,k+1:n) = C(k:m,k+1:n) - tao*v*(v'*C(k:m,k+1:n));
    
    % Update F (if exist) by the Householder transformation
        if p > 0
            F(k:m,:) = F(k:m,:) - tao*v*(v'*F(k:m,:));
        end
    end

    % Update the norm sof the column vectors
    for j = k+1:n
        s(j) = s(j) - C(k,j)^2;
    end
end

R = triu(C(1:n,1:n));


