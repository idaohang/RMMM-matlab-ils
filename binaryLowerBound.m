function [lb,sHat] = binaryLowerBound(R,y,s,k)
% lb = lowerBound(R,y,s) calculate a lower bound on the solution to the ILS
%      problem
%
% Input arguments:
%    R ---- n by n real nonsingular upper triangular matrix
%    y ---- n-dimensional real vector
%    s ---- n dimensional integer vector, the last k elements are a current
%           estimate
%    k ---- the number of elements that are fixed in s
% Output arguments:
%    lb ---- a lower bound on norm(R(1:k-1,1:k-1)s(1:k-1) - y(1:k-1))

    m = length(R);
    z = y(1:k-1) - R(1:k-1,k:m)*s(k:m)';
    
    %should have SVD as a parameter because it should be pre-computed to
    %avoid wasteful computations at each tree lvl
    [U S V] = svd(R(1:k-1,1:k-1));
    q = U'*z;
    r = rank(R(1:k-1,1:k-1));
    
    testSum = 0;
    for i=1:r
        testSum = testSum+(q(i)/S(i,i))^2;
    end
    
    if(testSum > k-0.25)
        % Use the secant method to solve for a positive value of lamda
        a = 0;
        b = 0.5;
        fa = 0;
        fb = 0;
        for i=1:r
            fa = fa + ((S(i,i)*q(i))/(S(i,i) + a))^2;
            fb = fb + ((S(i,i)*q(i))/(S(i,i) + b))^2;
        end
        fa = fa - k + 0.25;
        fb = fb - k + 0.25;

        %root will be stored in b after completion
        while(fb > 10^-8)
            if(a == b)
                break;
            end
            a = a - fa/((fb-fa)/(b-a));
            fa = 0;
            for i=1:r
                fa = fa + ((S(i,i)*q(i))/(S(i,i) + a))^2;
            end
            fa = fa - k + 0.25;

            if(a == b)
                break;
            end
            b = b - fb/((fa-fb)/(a-b));
            fb = 0;
            for i=1:r
                fb = fb + ((S(i,i)*q(i))/(S(i,i) + b))^2;
            end
            fb = fb - k + 0.25;           
        end
        
        %compute sHat, needed in next step...
        sHat = zeros(1,k-1);
        temp = zeros(1,k-1);
        for i=1:r
            sHat = sHat + ((S(i,i)*q(i))/(S(i,i) + b))*V(:,i)';
            temp = temp + ((b*q(i))/(S(i,i)^2 + b)) * V(:,i)';
        end
            
        eVals = eig(R(1:k-1,1:k-1)'*R(1:k-1,1:k-1));
        lb = (b + eVals(1)) * norm(round(sHat) - sHat)^2 + temp;

    else
        lb = 0;
        sHat = zeros(1,k-1);
        for i=1:r
            sHat(i) = (q(i)/S(i,i))*V(:,i);
        end
    end
end