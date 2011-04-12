runs = 500;
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
m = 20;
n = 20; 
%B = randn(m,n)./100 + eye(m,n);
B = randn(m,n);
z_true = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*10);
y = B*z_true + 0.6*randn(m,1);
p = 1;

[Q R] = qr(B);
y2 = Q'*y;

tic
[Zhat2,numExpanded2] = purebfsearchconstrained(R,y2,-10,10,100);
time2(i) = toc;
expand2(i) = numExpanded2;

tic
[Zhat,numExpanded] = searchconstrained(R,y2,-10,10,1);
time(i)=toc;
expand(i) = numExpanded;

%numExpanded
%numExpanded2

if(norm(Zhat - Zhat2')~=0)
    error('Wrong Answer!');
end

end