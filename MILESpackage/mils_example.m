
%
% A small example to run function mils.m
%

% Construct data
m = 7;
k = 2;
n = 3;
p = 3;

A = randn(m,k);
B = randn(m,n);
y = randn(m,1);

display('Three pairs of optimal least squares solutions')
[X,Z] = mils(A,B,y,p)