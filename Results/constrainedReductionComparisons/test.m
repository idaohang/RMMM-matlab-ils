runs = 1000000;
time3 = zeros(1,runs);
time1 = zeros(1,runs);
time2 = zeros(1,runs);
flops3 = zeros(1,runs);
flops1 = zeros(1,runs);
flops2 = zeros(1,runs);

sz = 30;
for i = 1:runs
    i
    %
    % A small example to run function ils.m
    %
    % Construct data
    m = sz;
    n = m; 
    lower = -10;
    upper = 10;
    %B = randn(m,n)./100 + eye(m,n);
    B = randn(m,n);
    %B = randILS(m,n,12);
    z_true = (-1*ones(n,1)).^(mod(round(rand(n,1)*10),2)+1).*round(rand(n,1)*upper);
    %y = (-1*ones(m,1)).^(mod(round(rand(m,1)*10),2)+1).*round(rand(m,1)*upper);
    y = B*z_true + 0.3*randn(m,1);
    p = 1;

    l = -1*(1:n);
    u = 1:n;

    %flops(0);
%     tic
    [P2 z2] = constrainedReduction(B,y,l,u);
%     time2(i) = toc;
    %flops2(i) = flops;

    %flops(0);
%     tic
    [P3 z3] = newConstrainedReduction(B,y,l,u);
%     time3(i) = toc; 
    %flops3(i) = flops;

    %flops(0);
%     tic
%     [Q R] = qr(B);
%     Q1 = Q(:,1:n);
%     y2 = Q1*Q1'*y ;
    [P z] = otherConstrainedReduction(B,y,l,u);
%     time1(i) = toc; 
%     flops1(i) = flops;

    if(norm(P2 - P) ~= 0)
         error('Wrong Answer!');
    end
    
    if(norm(P3-P) ~= 0)
         error('Wrong Answer!');
    end
end
