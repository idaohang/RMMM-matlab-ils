function rate = successRate(R,sigma)
    [m n] = size(R);
    rate = 1;
    for i = 1:n
        rate = rate * (2*normcdf(abs(R(i,i))/(2*sigma),0,1) - 1);
    end
end