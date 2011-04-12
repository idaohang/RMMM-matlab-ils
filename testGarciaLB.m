runs=500;
for i = 1:runs
    i
    %
    % A small example to run function ils.m
    %
    % Construct data
    m = 32;
    n = 32; 
    %B = randn(m,n)./100 + eye(m,n);
    B = randn(m,n);
    z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*10);
    y = B*z_true + 0.1*randn(m,1);
    p = 1;

    % Reduction
    [R,Z,y2] = reduction(B,y);
    % Search
    [Zhat,numExpanded] = search(R,y2,1);
    % Perform the unimodual transformation to obtain the optimal solutions
    Zhat = Z*Zhat;

    dist = norm(B*Zhat-y);

    realSoln = B\y;
    U = inv(B);
    for j = 1:n
        lb = abs(realSoln(j) - round(realSoln(j)))/norm(U(:,j));
        if(lb > dist)
            error('Lower Bound Incorrect!');
        end
    end

end
