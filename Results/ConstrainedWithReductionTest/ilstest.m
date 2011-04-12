runs = 100;
time = zeros(1,runs);
time4 = zeros(1,runs);
time6 = zeros(1,runs);
time7 = zeros(1,runs);
time8 = zeros(1,runs);
time10 = zeros(1,runs);
expand = zeros(1,runs);
expand4 = zeros(1,runs);
expand6 = zeros(1,runs);
expand7 = zeros(1,runs);
expand8 = zeros(1,runs);
expand10 = zeros(1,runs);
for i = 1:runs

%
% A small example to run function ils.m
%
% Construct data
m = 32;
n = 32; 
%B = randn(m,n)./100 + eye(m,n);
B = randn(m,n);
z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*10);
y = B*z_true + 0.5*randn(m,1);
p = 1;
i
lower = -10;
upper = 10;
% Reduction
% [R,Z,y2] = reduction(B,y);

% [Z R] = qr(B);
% y2 = Z'*y;
l = lower*ones(n,1);
u = upper*ones(n,1);
[R,y2,babai,l,u,P] = constrainedReduction(B,y,l,u);

% [Z R ind y2] = gradientReduction(B,y,z_true);

% Search

%SE Search
% tic
% [Zhat,numExpanded] = search(R,y2(1:n),p);
% time(i) = toc;
%SE Search with lb
% tic
% [Zhat2,numExpanded2] = searchlb(R,y2(1:n),p);
% toc
%SE Search with partial stopping radius
% tic
% [Zhat3,numExpanded3] = searchsr(R,y2(1:n),p);
% toc
%SE Search with single LB !MAY NOT GIVE CORRECT SOL'N!
% tic
% [Zhat3,numExpanded3] = searchsinglelb(R,y2(1:n),p);
% toc
% %BF Search
% tic
% [Zhat4, numExpanded4] = purebfsearch(R,y2(1:n),100000);
% time4(i) = toc;
% % %BF Search with lb
% tic
% [Zhat5, numExpanded5] = purebfsearchlb(R,y2(1:n),80000);
% toc

%BestFirst Box contrained Search
tic
[Zhat8,numExpanded8] = purebfsearchconstrained(R,y2(1:n),lower,upper,500000);
time8(i) = toc;
%A* Box contrained Search
% tic
% [Zhat9,numExpanded9] = purebfsearchconstrainedlb(R,y2(1:n),lower,upper,100000);
% toc
%SE Box constrained Search with LB
tic
[Zhat7,numExpanded7] = searchconstrainedlb(R,y2(1:n),lower,upper,p);
time7(i) = toc;
%SE Box contrained Search
tic
[Zhat6,numExpanded6] = searchconstrained(R,y2(1:n),lower,upper,p);
time6(i) = toc;


expand6(i) = numExpanded6;
expand7(i) = numExpanded7;
expand8(i) = numExpanded8;
% expand(i) = numExpanded;
% expand4(i) = numExpanded4;
%expand10(i) = numExpanded10;

% (numExpanded/t1)/(numExpanded2/t2)

% numExpanded
%numExpanded2
% numExpanded3
% numExpanded4
% numExpanded5
numExpanded6
numExpanded7
numExpanded8
% numExpanded9
%norm(Zhat-Zhat2)
% norm(Zhat - Zhat3)
% norm(Zhat - Zhat4')
% norm(Zhat - Zhat5')
%norm(Zhat - Zhat6)

if(norm(Zhat6 - Zhat7) > 0  || numExpanded6 + 1 < numExpanded7)
    error('Wrong Answer 7!');
end
if(norm(Zhat6 - Zhat8') > 0  || numExpanded6 + 1 < numExpanded8)
    error('Wrong Answer 8!');
end
% if(norm(Zhat - Zhat4') >0 || numExpanded + 1 < numExpanded4)
%     error('Wrong Answer 4!');
% end
% norm(Zhat6 - Zhat9')

% Perform the unimodual transformation to obtain the optimal solutions
% Zhat = Z*Zhat;
% Zhat2 = Z*Zhat2;
% Zhat3 = Z*Zhat3';
% Zhat4 = Z*Zhat4';
% norm(Zhat-Zhat2)
% norm(B*Zhat - y)
% %LowerBound
% eigVals = eig(R'*R);
% LinvB = inv(R) * y;
% lb = eigVals(1)*norm(round(LinvB) - LinvB)^2

%test
end
