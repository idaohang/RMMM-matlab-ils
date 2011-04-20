function [S,w]=HQR(A,b)
[m,n]=size(A);
r=min(m,n);
for j=1:r
   x=A(j:m,j);
   [u,tao]=house(x);
   A(j:m,j:n)=A(j:m,j:n)-tao*u*(u'*A(j:m,j:n));
   b(j:m)=b(j:m)-tao*u*(u'*b(j:m));
   if( j<m )
     A(j+1:m,j)=0;
  end
end
A(r+1:m,:)=0;
S=A(1:r,:);
w=b(1:r);
