
clear; clc;

n = 150;
maxIter = 30;
time = zeros(maxIter,3);
condNum = zeros(maxIter,3);
lb = -inf(n,1);
ub = inf(n,1);

for count = 1:maxIter
count
H = genMatrixAll(n,8);
x = round(rand(n,1)*10);
y = H*x + 0.7*randn(n,1);

[t_z1,z_ils] = matlabilstest(H,y);
time(count,1) = t_z1;

condH = cond(H);
condNum(count,1) = condH;
end

figure;
plot(time(1:maxIter,1), ':.r', 'LineWidth', 1);


