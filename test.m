runs = 100;
time3 = zeros(1,runs);
time1 = zeros(1,runs);
time2 = zeros(1,runs);
expand1 = zeros(1,runs);
expand2 = zeros(1,runs);
expand3 = zeros(1,runs);
checkSum = zeros(1,runs);
offDiagSum = zeros(1,runs);
rowSum = zeros(1,runs);
checkSum2 = zeros(1,runs);
offDiagSum2 = zeros(1,runs);
rowSum2 = zeros(1,runs);
babai1 = zeros(1,runs);
babai2 = zeros(1,runs);
babai3 = zeros(1,runs);
babaiCorrect1 = zeros(1,runs);
babaiCorrect2 = zeros(1,runs);
ilsCorrect1 = zeros(1,runs);
ilsCorrect2 = zeros(1,runs);
successRate1 = zeros(1,runs);
successRate2 = zeros(1,runs);
successRate3 = zeros(1,runs);
sdevs = zeros(1,runs);
snrs = zeros(1,runs);


sz = 40;
qam = 8;
m = sz;
n = m;
sigma2 = 0.9;
for i = 1:runs
    i
    %
    % A small example to run function ils.m
    %
    % Construct data
    B = randn(m,n);
    z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*ceil(rand(m,1)*sqrt(qam));
    lower = -10;
    upper = 10;
    sdevs(i) = sigma2;
    noiseVec = sigma2*(0 + randn(m,1));
    snrs(i) = norm(B*z_true)^2/(n*sigma2^2);
    snrs(i)
    y = B*z_true + noiseVec;
    
    l = -inf*(ones(1,n));
    u = inf*(ones(1,n));

    [R2 Z y2] = reduction(B,y);    
    
    tic;
    [zhat,numExpanded] = search(R2,y2,1);
    time1(i) = toc;
    numExpanded
    expand1(i) = numExpanded;
    successRate1(i) = successRate(R2,sigma2);
    
    [P z] = otherConstrainedReduction(R2,y2,l,u);
    [Q R3] = qr(R2(:,P));
    y3 = Q'*y2;
    tic;
    [zhat2,numExpanded2] = search(R3,y3,1);
    time2(i) = toc;
    numExpanded2
    expand2(i) = numExpanded2;
    successRate2(i) = successRate(R3,sigma2);
    [checkSum(i),rowSum(i),offDiagSum(i)] = checkLLL(R3);
    [sorted,idx] = sort(P);
    
    if(norm(zhat-zhat2(idx)) ~= 0)
        error('Wrong answer!');
    end
    
    [R4 Z4 y4] = testReduction(B,y,snrs(i)^(0.52),sigma2);
    tic;
    [zhat3,numExpanded3] = search(R4,y4,1);
    time3(i) = toc;
    numExpanded3
    expand3(i) = numExpanded3;
    successRate3(i) = successRate(R4,sigma2);
    [checkSum2(i),rowSum2(i),offDiagSum2(i)] = checkLLL(R4);
    
    if(norm(Z*zhat-Z4*zhat3) ~= 0)
        error('Test Reduction Wrong Answer!');
    end
    
    babai1(i) = norm(B*Z*babai(R2,y2) - y);
    if(norm(Z*babai(R2,y2) - z_true) == 0)
        babaiCorrect1(i) = 1;
    end
    
    babai2temp = babai(R3,y3);
    babai2(i) = norm(B*Z*babai2temp(idx)-y);
    if(norm(Z*babai2temp(idx) - z_true) == 0)
        babaiCorrect2(i) = 1;
    end
    
end

