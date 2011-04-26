runs = 100;
checkSum = zeros(6,9,runs);
offDiagSum = zeros(6,9,runs);
rowSum = zeros(6,9,runs);
numSwapped = zeros(6,9,runs);
avgAngle = zeros(6,9,runs);

sz = 40;
qam = 8;
sigmas = 0.1:0.1:0.9;
m = sz;
n = m;
for k = 20:5:45
    m = k;
    n = k;
    for j = 1:9
        for i = 1:runs
            fprintf('%i %i %i\n',k,j,i);
            sigma = sigmas(j);
            %
            % A small example to run function ils.m
            %
            % Construct data
            lower = -10;
            upper = 10;
            B = randn(m,n);
            z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*ceil(rand(m,1)*sqrt(qam));
            %y = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*upper);
            y = B*z_true + sigma*randn(m,1);
            %y = randn(m,1)*20;
            p = 1;

            l = -inf*(ones(1,n));
            u = inf*(ones(1,n));

            [R2 Z y2] = reduction(B,y);

            [P z] = otherConstrainedReduction(R2,y2,l,u);
            [Q R3] = qr(R2(:,P));
            y3 = Q'*y2;

            [checkSum(k,j,i),rowSum(k,j,i),offDiagSum(k,j,i),avgAngle(k,j,i)] = checkLLL(R3);

            [sorted,idx] = sort(P);
            numSwapped(k,j,i) = size(find(idx-P),2);


        end
    end
end

