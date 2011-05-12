runs = 25;
time3 = zeros(3,runs);
time1 = zeros(3,runs);
time2 = zeros(3,runs);
expand1 = zeros(3,runs);
expand2 = zeros(3,runs);
expand3 = zeros(3,runs);
checkSum = zeros(3,runs);
offDiagSum = zeros(3,runs);
rowSum = zeros(3,runs);
checkSum2 = zeros(3,runs);
offDiagSum2 = zeros(3,runs);
rowSum2 = zeros(3,runs);
babai1 = zeros(3,runs);
babai2 = zeros(3,runs);
babai3 = zeros(3,runs);
babaiCorrect1 = zeros(3,runs);
babaiCorrect2 = zeros(3,runs);
ilsCorrect1 = zeros(3,runs);
ilsCorrect2 = zeros(3,runs);
successRate1 = zeros(3,runs);
successRate2 = zeros(3,runs);
successRate3 = zeros(3,runs);
sdevs = zeros(3,runs);
snrs = zeros(3,runs);


sz = 40;
qam = 8;
m = sz;
n = m;
B = randn(m,n)/10;
z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*ceil(rand(m,1)*sqrt(qam));
tempSNR = [5,8,10];
for k=1:3
for i = 1:runs
    i
    %
    % A small example to run function ils.m
    %
    % Construct data
    lower = -10;
    upper = 10;
    
    sigma2 = sqrt((norm(B*z_true)^2)/tempSNR(k))/n
    sdevs(i) = sigma2;
    tempSNR(k)
    %snrs(i) = norm(B*z_true)^2/(n*sigma2)^2;
    y = B*z_true + sigma2*randn(m,1);
    
    l = -inf*(ones(1,n));
    u = inf*(ones(1,n));

    [R2 Z y2] = reduction(B,y);    
    
    tic;
    [zhat,numExpanded] = search(R2,y2,1);
    time1(k,i) = toc;
    numExpanded
    expand1(k,i) = numExpanded;
    successRate1(k,i) = successRate(R2,sigma2);
    
    [P z] = otherConstrainedReduction(R2,y2,l,u);
    [Q R3] = qr(R2(:,P));
    y3 = Q'*y2;
    tic;
    [zhat2,numExpanded2] = search(R3,y3,1);
    time2(k,i) = toc;
    numExpanded2
    expand2(k,i) = numExpanded2;
    successRate2(k,i) = successRate(R3,sigma2);
    [checkSum(k,i),rowSum(k,i),offDiagSum(k,i)] = checkLLL(R3);
    [sorted,idx] = sort(P);
    
    if(norm(zhat-zhat2(idx)) ~= 0)
        error('Wrong answer!');
    end
    
    [R4 Z4 y4] = testReduction(B,y,(1/sigma2)^3*3);
    tic;
    [zhat3,numExpanded3] = search(R4,y4,1);
    time3(k,i) = toc;
    expand3(k,i) = numExpanded3;
    successRate3(k,i) = successRate(R4,sigma2);
    [checkSum2(k,i),rowSum2(k,i),offDiagSum2(k,i)] = checkLLL(R4);
    
    if(norm(Z*zhat-Z4*zhat3) ~= 0)
        error('Test Reduction Wrong Answer!');
    end
    
    babai1(i) = norm(B*Z*babai(R2,y2) - y);
    if(norm(Z*babai(R2,y2) - z_true) == 0)
        babaiCorrect1(k,i) = 1;
    end
    
    babai2temp = babai(R3,y3);
    babai2(k,i) = norm(B*Z*babai2temp(idx)-y);
    if(norm(Z*babai2temp(idx) - z_true) == 0)
        babaiCorrect2(k,i) = 1;
    end
    
end
end

