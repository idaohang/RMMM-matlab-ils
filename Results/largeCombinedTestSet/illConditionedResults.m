sizes = 35;
noises = [0.1;0.3];
runs=100;

babaiDists = zeros(1,2,runs);
volumes = zeros(1,2,runs);
volumeRatios=zeros(1,2,runs);

time = zeros(1,2,runs,25);
expand = zeros(1,2,runs,25);

for i = 1:1
    for j = 1:2
        for p = 1:runs
            % Construct data
            m = sizes(i);
            n = sizes(i);
            %B = randn(m,n)./100 + eye(m,n);
            B = genMatrixAll(sizes,10);
            z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*10);
            y = B*z_true + noises(j)*randn(m,1);

            % [Q R] = qr(B);
            % y2 = Q'*y;
            [R Z y2] = reduction(B,y);
            babaiDists(i,j,p) = norm(R*babai(R,y2)-y2)^2;
            volumes(i,j,p) = abs(det(R));
            if(mod(n,2) == 0)
                sphereVolume = ((pi^(n/2))/factorial(n/2))*(babaiDists(i,j,p)^n);
                volumeRatios(i,j,p) = sphereVolume/volumes(i,j,p);
            else
                sphereVolume = ((2^((n+1)/2))*(pi^((n-1)/2)))/factorial(factorial(n));
                volumeRatios(i,j,p) = sphereVolume/volumes(i,j,p);
            end
            
            Zhat = zeros(n,1);
            tempCheck = zeros(n,1);
            for k = 2:2:sizes(i)
                fprintf('%i %i %i %i\n',i,j,p,k/2);
                if(k~=4 && norm(tempCheck-Zhat) ~= 0)
                    error('Different answer with changing k!');
                end
                tempCheck = Zhat;
                tic
                [Zhat,numExpanded] = bfsdfscombined(R,y2,k,10,100,10000000);
                time(i,j,p,k/2) = toc;
                expand(i,j,p,k/2) = numExpanded;
            end
        end
    end
end
