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
% Version 1.1,  January 2011.


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
s = zeros(1,n);
piv = zeros(n,1);

for k = 1:n

    % Compute the 2-norm of each column of C(k:m,k:n)
    for j = k:n
        s(j) = norm(C(k:m,j));  
    end
    
    % Find the column with minimum 2-norm 
    [minnorm, i] = min(s(k:n));
    q = i+k-1;
    piv(k) = q;
    
    % Column interchange
    if q > k
        C(:,[k,q]) = C(:,[q,k]);
    end

    % Compute and apply the Householder transformation  I-tau*v*v'
    if norm(C(k+1:m,k)) > 0, % otherwise no Householder transformation is needed
	    v = C(k:m,k);
        rho = minnorm;     
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
end

R = triu(C(1:n,1:n));


