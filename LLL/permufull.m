function [R,y,y2]=permufull(R,y,y2,j,k)
%Swap columns i and k where k > j, then bring R back to upper triangular with a
%product of givens rotations
if(k<j)
    temp = k;
    k=j;
    j=temp;
end

m = size(R,2);
n=k;
Rj = R(:,j);
R(:,j) = R(:,n);
R(:,n) = Rj;

for i=n-1:-1:j
    r=sqrt(R(i,j)*R(i,j)+R(i+1,j)*R(i+1,j));
    c=R(i,j)/r;
    s=R(i+1,j)/r;
    G=[c,s; -s,c];
    R([i,i+1],i:m)=G*R([i,i+1],i:m);
    R(i,j)=r;
    R(i+1,j)=0;       
    y([i,i+1]) = G*y([i,i+1]);
    y2([i,i+1]) = G*y2([i,i+1]);
    
end

for i = j+1:m-1
    r=sqrt(R(i,i)*R(i,i)+R(i+1,i)*R(i+1,i));
    c=R(i,i)/r;
    s=R(i+1,i)/r;
    G=[c,s; -s,c];
    R(i,i)=r;
    R(i+1,i)=0;
    R([i,i+1],i+1:m)=G*R([i,i+1],i+1:m);
    y([i,i+1]) = G*y([i,i+1]);
    y2([i,i+1]) = G*y2([i,i+1]);

end

