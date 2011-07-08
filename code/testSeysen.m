n = 70;
sigma = 0.5;

maxIter = 100;
time1 = zeros(1,maxIter);
time2 = zeros(1,maxIter);
condNum = zeros(1,maxIter);
lb = -inf(n,1);
ub = inf(n,1);
babBetter = zeros(1,maxIter);

for count = 1:maxIter
count
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
babai1(count) = norm(H*Z2*babai(R2,y2) - y);
if(norm(Z2*babai(R2,y2) - x) == 0)
    babaiCorrect1(count) = 1;
end

babai2temp = babai(R3,y3);
babai2(count) = norm(H*Z2*babai2temp(idx)-y);
if(norm(Z2*babai2temp(idx) - x) == 0)
    babaiCorrect2(count) = 1;
end
babai1(count)
babai2(count)


tic
[z_ils] = mexsearch(R2,y2);
t_z1 = toc
z_ils = Z2*z_ils;


if(babai2(count) < babai1(count))
    babBetter(count) = 1;
    tic;
    [z_ils2] = mexsearch(R3,y3);
    t_z2 = toc
    [~, idx] = sort(P);
    z_ils2 = Z2*z_ils2(idx);
else
    tic;
    [z_ils2] = mexsearch(R2,y2);
    t_z2 = toc 
    z_ils2 = Z2*z_ils2;
end


%z_ils2 = Z3*z_ils2;


if(norm(z_ils-z_ils2) ~= 0)
    error('Wrong Answer!');
end

time1(count) = t_z1;
time2(count) = t_z2;

end

mean(time1)
mean(time2)

