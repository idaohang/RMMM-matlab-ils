sizes = 10:5:50;
noises = [0.05;0.1;0.3;0.4;0.5;0.7;0.8];
runs=200;

time1 = zeros(9,7,runs);
expand1 = zeros(9,7,runs);

time2 = zeros(9,7,runs);
expand2 = zeros(9,7,runs);

timePermute = zeros(9,7,runs);

for i = 1:size(sizes,2)
    for j = 1:size(noises,1)
        for p = 1:runs
            fprintf('%i %i %i\n',i,j,p);
            % Construct data
            m = sizes(i);
            n = sizes(i);
            %B = randn(m,n)./100 + eye(m,n);
            B = randn(m,n);
            z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*10000);
            y = B*z_true + noises(j)*randn(m,1);
            l = -inf*(ones(1,n));
            u = inf*(ones(1,n));

            [R Z y2] = reduction(B,y);
            
            tic;
            [zhat,numExpanded] = search(R,y2,1);
            time1(i,j,p) = toc;
            
            tic;
            [P z] = otherConstrainedReduction(R,y2,l,u);
            [Q R2] = qr(R(:,P));
            y3 = Q'*y2;
            timePermute(i,j,p) = toc;

            tic;
            [zhat2,numExpanded2] = search(R2,y3,1);
            time2(i,j,p) = toc;
            
            [sorted,idx] = sort(P);
            if(norm(zhat-zhat2(idx)) ~= 0)
                error('Wrong answer!');
            end
            
            expand1(i,j,p) = numExpanded;
            expand2(i,j,p) = numExpanded2;
        end
    end
end