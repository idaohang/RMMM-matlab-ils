function [R,Z,y,y2]=swap(R,k,Z,y,y2)

[m,n]=size(R);
Rk=R(:,k);
R(:,k)=R(:,k+1);
R(:,k+1)=Rk;

Z(:,[k,k+1])=Z(:,[k+1,k]);

r=sqrt(R(k,k)^2+R(k+1,k)^2);

c=R(k,k)/r;
s=R(k+1,k)/r;
G=[c,s; -s,c];

R(k,k)=r;
R(k+1,k)=0;
R([k,k+1],k+1:n)=G*R([k,k+1],k+1:n);

y([k,k+1]) = G*y([k,k+1]);
y2([k,k+1]) = G*y2([k,k+1]);
