runs = 1000;
time3 = zeros(2,runs);
time1 = zeros(2,runs);
time2 = zeros(2,runs);
time4 = zeros(2,runs);
expand1 = zeros(2,runs);
expand2 = zeros(2,runs);
expand3 = zeros(2,runs);
expand4 = zeros(2,runs);
checkSum = zeros(2,runs);
offDiagSum = zeros(2,runs);
rowSum = zeros(2,runs);
checkSum2 = zeros(2,runs);
offDiagSum2 = zeros(2,runs);
rowSum2 = zeros(2,runs);
checkSum3 = zeros(2,runs);
offDiagSum3 = zeros(2,runs);
rowSum3 = zeros(2,runs);
babai1 = zeros(2,runs);
babai2 = zeros(2,runs);
babai3 = zeros(2,runs);
babaiCorrect1 = zeros(2,runs);
babaiCorrect2 = zeros(2,runs);
ilsCorrect1 = zeros(2,runs);
ilsCorrect2 = zeros(2,runs);
successRate1 = zeros(2,runs);
successRate2 = zeros(2,runs);
successRate3 = zeros(2,runs);
successRate4 = zeros(2,runs);
numSwaps = zeros(2,runs);


sz = 40;
qam = 8;
sigmas = [0.4,0.9];
m = sz;
n = m;
for nz = 1:2
    sigma = sigmas(nz);
    for i = 1:runs
        i
        %
        % A small example to run function ils.m
        %
        % Construct data
        lower = -10;
        upper = 10;
        B = randn(m,n);
        %B = randILS(m,n,6);
        %B = randILS(m,n,12);
        z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*ceil(rand(m,1)*sqrt(qam));
        %y = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*upper);
        y = B*z_true + sigma*randn(m,1);
        %y = randn(m,1)*20;
        p = 1;

        l = -inf*(ones(1,n));
        u = inf*(ones(1,n));

        [R2 Z y2] = reduction(B,y);


        tic;
        [zhat,numExpanded] = search(R2,y2,1);
        time1(nz,i) = toc;
        expand1(nz,i) = numExpanded;
        successRate1(nz,i) = successRate(R2,sigma);

        [P z] = otherConstrainedReduction(R2,y2,l,u);
        [Q R3] = qr(R2(:,P));
        y3 = Q'*y2;
        tic;
        [zhat2,numExpanded2] = search(R3,y3,1);
        time2(nz,i) = toc;
        expand2(nz,i) = numExpanded2;
        successRate2(nz,i) = successRate(R3,sigma);
        [checkSum(nz,i),rowSum(nz,i),offDiagSum(nz,i)] = checkLLL(R3);
        [sorted,idx] = sort(P);

        if(norm(zhat-zhat2(idx)) ~= 0)
            error('Wrong answer!');
        end

        numSwaps(nz,i) = sum(P~=1:n);
    %     [R4 Z4 y4] = testReduction(B,y,4);
    %     tic;
    %     [zhat3,numExpanded3] = search(R4,y4,1);
    %     time3(i) = toc;
    %     expand3(i) = numExpanded3;
    %     successRate3(i) = successRate(R4,sigma);
    %     [checkSum2(i),rowSum2(i),offDiagSum2(i)] = checkLLL(R4);
    %     
    %     if(norm(Z*zhat-Z4*zhat3) ~= 0)
    %         error('Test Reduction Wrong Answer!');
    %     end
    %     
    %     [R5 Z5 y5] = testReduction2(B,y,4);
    %     tic;
    %     [zhat4,numExpanded4] = search(R5,y5,1);
    %     time4(i) = toc;
    %     expand4(i) = numExpanded4;
    %     successRate4(i) = successRate(R5,sigma);
    %     [checkSum3(i),rowSum3(i),offDiagSum3(i)] = checkLLL(R5);
    %     
    %     if(norm(Z*zhat-Z5*zhat4) ~= 0)
    %         error('Test Reduction2 Wrong Answer!');
    %     end
    %     
    %     babai1(i) = norm(B*Z*babai(R2,y2) - y);
    %     if(norm(Z*babai(R2,y2) - z_true) == 0)
    %         babaiCorrect1(i) = 1;
    %     end
    %     
    %     babai2temp = babai(R3,y3);
    %     babai2(i) = norm(B*Z*babai2temp(idx)-y);
    %     if(norm(Z*babai2temp(idx) - z_true) == 0)
    %         babaiCorrect2(i) = 1;
    %     end
    %     expand1(i) = numExpanded;
    %     expand2(i) = numExpanded2;

    end

end
