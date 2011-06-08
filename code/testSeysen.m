
clear; clc;

n = 15;
maxIter = 100;
time = zeros(maxIter,3);
condNum = zeros(maxIter,3);
lb = -inf(n,1);
ub = inf(n,1);

for count = 1:maxIter
count
H = genMatrixAll(n,10);
x = round(rand(n,1)*10);
y = H*x + 0.01*randn(n,1);

[t_z1,z_ils] = matlabilstest(H,y);
[t_z2,z_ils2] = meschachmatlab(H,y);

if(norm(z_ils-z_ils2) ~= 0)
    error('Wrong Answer!');
end

time(count,1) = t_z1;
time(count,2) = t_z2;
condH = cond(H);
condNum(count,1) = condH;
end

figure;
plot(1:maxIter,time(1:maxIter,1), '--.r',1:maxIter,time(1:maxIter,2), 'xb-');

