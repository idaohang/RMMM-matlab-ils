sizes = 10:5:50;
noises = [0.05;0.1;0.3;0.4;0.5;0.7;0.8];
runs=200;

time1 = zeros(9,7,runs);
expand1 = zeros(9,7,runs);
babai1 = zeros(9,7,runs);
babaiCorrect1 = zeros(9,7,runs);

time2 = zeros(9,7,runs);
expand2 = zeros(9,7,runs);
babai2 = zeros(9,7,runs);
babaiCorrect2 = zeros(9,7,runs);

time3 = zeros(9,7,runs);

babaiBetter = zeros(9,7,runs);

checkSum = zeros(9,7,runs);
offDiagSum = zeros(9,7,runs);
rowSum = zeros(9,7,runs);

timePermute = zeros(9,7,runs);
qam=64;
for i = 1:size(sizes,2)
    for j = 1:size(noises,1)
        for p = 1:runs
            fprintf('%i %i %i\n',i,j,p);
            % Construct data
            m = sizes(i);
            n = sizes(i);
            %B = randn(m,n)./100 + eye(m,n);
            B = randn(m,n);
            z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*ceil(rand(m,1)*sqrt(qam));
            y = B*z_true + noises(j)*randn(m,1);
            l = -inf*(ones(1,n));
            u = inf*(ones(1,n));

            [R Z y2] = reduction(B,y);
                        
            babai1(i,j,p) = norm(B*Z*babai(R,y2) - y);
            if(norm(Z*babai(R,y2) - z_true) == 0)
                babaiCorrect1(i,j,p) = 1;
            end
            
            tic;
            [zhat,numExpanded] = search(R,y2,1);
            time1(i,j,p) = toc;
            
            tic;
            [P z] = otherConstrainedReduction(R,y2,l,u);
            [Q R2] = qr(R(:,P));
            y3 = Q'*y2;
            timePermute(i,j,p) = toc;
            [checkSum(i,j,p),rowSum(i,j,p),offDiagSum(i,j,p)] = checkLLL(R2);   

            [sorted,idx] = sort(P);
            babai2temp = babai(R2,y3);
            babai2(i,j,p) = norm(B*Z*babai2temp(idx)-y);
            if(norm(Z*babai2temp(idx) - z_true) == 0)
                babaiCorrect2(i) = 1;
            end
            
            if(babai2(i,j,p) < babai1(i,j,p))
                babaiBetter(i,j,p) = 1;
                tic;
                [zhat2,numExpanded2] = search(R2,y3,1);
                time2(i,j,p) = toc;
                
                tic;
                [zhat3,numExpanded3] = searchExtra(R,y2,babai2(i,j,p)^2+10^-5,0);
                time3(i,j,p) = toc;
            else
                idx = 1:n;
                zhat2 = zhat;
                zhat3 = zhat;
                time2(i,j,p) = time1(i,j,p);
                numExpanded2 = numExpanded;
                
                time3(i,j,p) = time1(i,j,p);
                numExpanded3 = numExpanded;
            end
            
            if(norm(zhat-zhat2(idx)) ~= 0)
                error('Wrong answer!');
            end
            if(norm(zhat-zhat3)~=0)
                error('Wrong Answer!');
            end
            
            expand1(i,j,p) = numExpanded;
            expand2(i,j,p) = numExpanded2;
        end
    end
end