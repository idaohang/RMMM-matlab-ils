function [Zhat,numExpanded] = searchlb(R,y,p)
%
% Zhat = search(R,y,p) produces p optimal solutions to the upper triangular
%        integer least squares problem min_{z}||y-Rz|| by a search algorithm.
%
% Input arguments:
%    R ---- n by n real nonsingular upper triangular matrix
%    y ---- n-dimensional real vector
%    p ---- the number of optimal solutions and its default value is 1
%
% Output arguments:
%    Zhat - n by p integer matrix (in double precision). Its j-th column 
%           s the j-th optimal solution, i.e., its residual norm is the j-th
%           smallest, so ||y-R*Zhat(:,1)|| <= ...<= ||y-R*Zhat(:,p)||

% Copyright (c) 2006. Xiao-Wen Chang and Tianyang Zhou
% Version 1.0, October 2006.


% Check input arguments
if nargin < 2 % input error
    error('Not enough input arguments!')
end

if nargin < 3
    p = 1;
end

if p <= 0 % input error
    error('Third input argument must be an integer bigger than 0!');
end

[m,n] = size(R);
[n2,n3] = size(y);
if m ~= n | n ~= n2 | n3 ~= 1  % input error
    error('Input arguments have a matrix dimension error!')
end

% Initialization
z = zeros(n,1); % the current point
c = zeros(n,1); % c(k)=(y(k)-R(k,k+1:n)z(k+1:n))/R(k,k)
d = zeros(n,1); % d(k) is used to compute z(k) 
prsd = zeros(n,1); % partial squared residual norm for z
                   % prsd(k)=norm(y(k+1:n)-R(k+1:n,k+1:n)z(k+1:n))^2
S = zeros(n,n+1); % S(:,k) = R(:,k:n)*z(k:n), k=1:n
Zhat = zeros(n,p); % the p candidate solutions (or points) 
rsd = zeros(p,1); % squared residual norms of the p candidate solutions

beta = inf;  % the initial ellipsoid bound
ncand = 0;   % the initial number of candidate solutions

%Preprocessing for lower bound
[eVals F] = preProcessLB(R);
f = zeros(m,m);

c(n) = y(n)/R(n,n);
z(n) = round(c(n));
gamma = R(n,n)*(c(n)-z(n));
if c(n) > z(n)
    d(n) = 1;
else
    d(n) = -1;
end

k = n;
numExpanded=0;

while 1
    newprsd = prsd(k) + gamma*gamma;
    %Compute a lowerbound on the subproblem on the next level down
    lb = 0;
    if(newprsd < beta)
        if(k ~= 1 && k~= n)  
            %OLD SLOW WAY
            %f(1:k-1,k-1) = f(1:k-1,k) + (R(k,k+1:n)*z(k+1:m))*F(1:k-1,k) - F(1:k-1,k)*y(k) - FR(1:k-1,k)*z(k);
            %NEW WAY
            f(1:k-1,k-1) = f(1:k-1,k) - F(1:k-1,k)*gamma;
            lb = eVals(k-1)*norm((round(f(1:k-1,k-1)) - f(1:k-1,k-1)))^2;
        else
            if(k == n)
                f(1:m-1,m-1) = F(1:m-1,1:m-1)*(y(1:m-1) - R(1:m-1,m)*z(m));
                lb = eVals(m-1)*norm((round(f(1:m-1,m-1)) - f(1:m-1,m-1)))^2;
            end
        end
    end
    
    if newprsd < (beta - lb)
        if k ~= 1 % move to level k-1
            numExpanded=numExpanded+1;
            S(1:k,k) = R(1:k,k)*z(k) + S(1:k,k+1);
            k = k - 1;
            prsd(k) = newprsd;
            c(k) = (y(k)-S(k,k+1))/R(k,k);
            z(k) = round(c(k));
            gamma = R(k,k)*(c(k)-z(k));
            if c(k) > z(k) 
                d(k) = 1;
            else
                d(k) = -1;
            end
        else % a new point is found, update the set of candidate solutions
            if ncand < p   % add the new point
                numExpanded=numExpanded+1;
                ncand = ncand + 1;
                Zhat(:,ncand) = z;
                rsd(ncand) = newprsd;
                if ncand == p
                    beta = rsd(p);
                end
            else   % insert the new point and remove the worst one
                numExpanded=numExpanded+1;
                i = 1;
                while i < p && rsd(i) <= newprsd
                    i = i + 1;
                end
                Zhat(:,i:p) = [z, Zhat(:,i:p-1)];
                rsd(i:p) = [newprsd; rsd(i:p-1)];
                beta = rsd(p);
            end

            z(1) = z(1) + d(1);
            gamma = R(1,1)*(c(1)-z(1));
            if d(1) > 0
                d(1) = -d(1) - 1;
            else
                d(1) = -d(1) + 1;
            end
        end
    else  
        if (k == n)   % not sure about this line... used to be "if (k == n && sigma(k) >= beta)"..no idea where sigma(k) came from...
            break
        else  % move back to level k+1
            if(newprsd > beta)
                k = k + 1;
            end
            z(k) = z(k) + d(k);
            gamma = R(k,k)*(c(k)-z(k));
            if d(k) > 0
                d(k) = -d(k) - 1;
            else
                d(k) = -d(k) + 1;
            end
        end 
    end
end

