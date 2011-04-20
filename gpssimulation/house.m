function [u,tao] = house(x)
n=length(x);
norm2x=x'*x;
si=-sign(x(1));
p=si*sqrt(norm2x);
u=[x(1)-p;x(2:n)];
tao=-1/(p*(x(1)-p));
