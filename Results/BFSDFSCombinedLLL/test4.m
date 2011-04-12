runs = 100;
time3 = zeros(1,runs);
expand3 = zeros(1,runs);
time = zeros(1,runs);
expand = zeros(1,runs);
expand2 = zeros(1,runs);
time2 = zeros(1,runs);

for i = 1:runs
i
%
% A small example to run function ils.m
%
% Construct data
m = 45;
n = 45; 
%B = randn(m,n)./100 + eye(m,n);
B = randn(m,n);
z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*10);
y = B*z_true + 0.8*randn(m,1);
p = 1;

% [Q R] = qr(B);
% y2 = Q'*y;
[R Z y2] = reduction(B,y);


tic
[Zhat3,numExpanded3] = purebfsearch(R,y2,3000000);
time3(i) = toc;
expand3(i) = numExpanded3;

tic
[Zhat2,numExpanded2] = bfsdfscombined(R,y2,32,10,100,3000000);
time2(i) = toc;
expand2(i) = numExpanded2;

tic
[Zhat,numExpanded] = search(R,y2,1);
time(i)=toc;
expand(i) = numExpanded;

time(i)
time2(i)
time3(i)
numExpanded
numExpanded2
numExpanded3

if(norm(Zhat-Zhat2)~=0)
    error('Wrong Answer!');
end

end