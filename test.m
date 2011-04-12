runs = 10000;
time3 = zeros(1,runs);
time1 = zeros(1,runs);
time2 = zeros(1,runs);
expand1 = zeros(1,runs);
expand2 = zeros(1,runs);
expand3 = zeros(1,runs);
checkSum = zeros(1,runs);
offDiagSum = zeros(1,runs);
rowSum = zeros(1,runs);
babai1 = zeros(1,runs);
babai2 = zeros(1,runs);

sz = 3;
qam = 64;
for i = 1:runs
    i
    %
    % A small example to run function ils.m
    %
    % Construct data
    m = sz;
    n = m; 
    lower = -10;
    upper = 10;
    %B = randn(m,n);
    B = randILS(m,n,6);
    %B = randILS(m,n,12);
    z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*ceil(rand(m,1)*sqrt(qam));
    %y = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*upper);
    %y = B*z_true + 0.5*randn(m,1);
    y = randn(m,1)*20;
    p = 1;
    
    l = -inf*(ones(1,n));
    u = inf*(ones(1,n));

    
    [R2 Z y2] = reduction(B,y);
    
    tic;
    [zhat,numExpanded] = search(R2,y2,1);
    time1(i) = toc;
    
    [P z] = otherConstrainedReduction(R2,y2,l,u);
    [Q R3] = qr(R2(:,P));
    y3 = Q'*y2;
    tic;
    [zhat2,numExpanded2] = search(R3,y3,1);
    time2(i) = toc;
    
    [checkSum(i),rowSum(i),offDiagSum(i)] = checkLLL(R3);

    [sorted,idx] = sort(P);
    
    if(norm(zhat-zhat2(idx)) ~= 0)
        error('Wrong answer!');
    end
    
    babai1(i) = norm(B*Z*babai(R2,y2) - y);
    babai2temp = babai(R3,y3);
    babai2(i) = norm(B*Z*babai2temp(idx)-y);

    expand1(i) = numExpanded;
    expand2(i) = numExpanded2;
    
    if(numExpanded > numExpanded2 + 3 && numExpanded2 > 3)
        error('test');
    end

end
