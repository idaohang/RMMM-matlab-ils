function [R,y,L]=permuTest(R,y,j,L)
%Swap columns n and j, then bring R back to upper triangular with the
%product of givens rotations G. Also apply the givens rotations to a vector
%y and optionally a lower triangular matrix L.
[m,n] = size(R);

tempCol = R(:,j);
for i=j:n-1
    R(:,i) = R(:,i+1);
end
R(:,n) = tempCol;

if(nargin == 4)
    tempCol = L(:,j);
    for i=j:n-1
        L(:,i) = L(:,i+1);
    end
    L(:,n) = tempCol;
end

u=j;
for i=u:n-1
    j=i;
    r=sqrt(R(i,j)*R(i,j)+R(i+1,j)*R(i+1,j));
    c=R(i,j)/r;
    s=R(i+1,j)/r;
    G=[c,s; -s,c];
    R([i,i+1],i:n)=G*R([i,i+1],i:n);
    %addflops(flops_mul(G,R([i,i+1],i:n)));

    R(i,j)=r;
    R(i+1,j)=0;
       
    y([i,i+1]) = G*y([i,i+1]);
    %addflops(flops_mul(G,y([i,i+1])));
    
    if(nargin==4)
        L([i,i+1],[1:i,n]) = G*L([i,i+1],[1:i,n]);
        %addflops(flops_mul(G,L([i,i+1],[1:i,n])));
    end
end
