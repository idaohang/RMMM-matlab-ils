n = 45;
sigmas = [0.6];

maxIter = 1000;
time1 = zeros(length(sigmas),maxIter);
time2 = zeros(length(sigmas),maxIter);
time3 = zeros(length(sigmas),maxIter);
time4 = zeros(length(sigmas),maxIter);
condNum = zeros(length(sigmas),maxIter);
lb = -inf(n,1);
ub = inf(n,1);
babBetter = zeros(length(sigmas),maxIter);

outerCount=0;
for sigma = sigmas
    outerCount=outerCount+1;
    for count = 1:maxIter
    outerCount,count
    H = genMatrixAll(n,11);
    x = (-1*ones(n,1)).^(mod(round(rand(n,1)*10),2)+1).*round(rand(n,1)*10);
    snr = norm(H*x)^2/(n*sigma^2);
    v = sigma*randn(n,1);
    y = H*x + v;

    l = -inf*ones(n,1);
    u = inf*ones(n,1);
    [R2 Z2 y2] = reduction(H,y);
    P = otherConstrainedReduction(R2,y2,lb,ub);
    [Q3 R3] = qr(R2(:,P));
    y3 = Q3'*y2;
    %[R3 Z3 y3] = testReduction(H,y,snr^(0.5));

    [sorted,idx] = sort(P);
    babai1(outerCount,count) = norm(H*Z2*babai(R2,y2) - y);
    if(norm(Z2*babai(R2,y2) - x) == 0)
        babaiCorrect1(outerCount,count) = 1;
    end

    babai2temp = babai(R3,y3);
    babai2(outerCount,count) = norm(H*Z2*babai2temp(idx)-y);
    if(norm(Z2*babai2temp(idx) - x) == 0)
        babaiCorrect2(outerCount,count) = 1;
    end
    babai1(outerCount,count)
    babai2(outerCount,count)


    tic;
    [z_ils] = searchExtra(R2,y2,babai1(outerCount,count)^2+10^-6,0);
    t_z1 = toc
    z_ils = Z2*z_ils;

    tic;
        [z_ils4] = searchExtra(R3,y3,babai2(outerCount,count)^2+10^-6,0);
    t_z4 = toc
    [~, idx] = sort(P);
    z_ils4 = Z2*z_ils4(idx);

    if(babai2(outerCount,count) < babai1(outerCount,count))
        babBetter(outerCount,count) = 1;
        tic;
        [z_ils2] = searchExtra(R2,y2,babai2(outerCount,count)^2+10^-6,0);
        t_z2 = toc

        tic;
            [z_ils3] = searchExtra(R3,y3,babai2(outerCount,count)^2+10^-6,0);
        t_z3 = toc
        z_ils3 = Z2*z_ils3(idx);

        z_ils2 = Z2*z_ils2;
    else
        z_ils2 = z_ils;
        t_z2 = t_z1;
        %z_ils2 = Z2*z_ils2;
        z_ils3 = z_ils2;
        t_z3 = t_z2;
    end


    %z_ils2 = Z3*z_ils2;


    if(norm(z_ils-z_ils2) ~= 0 || norm(z_ils - z_ils3)~= 0 || norm(z_ils4 - z_ils) ~= 0 )
        error('Wrong Answer!');
    end

    
    if(t_z1 - 10 > t_z2)
        error('Got one');
    end
    time1(outerCount,count) = t_z1;
    time2(outerCount,count) = t_z2;
    time3(outerCount,count) = t_z3;
    time4(outerCount,count) = t_z4;

    
    end

end

