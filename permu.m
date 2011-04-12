function [R,y,L]=permu(R,y,j,L)
%Swap columns n and j, then bring R back to upper triangular with the
%product of givens rotations G. Also apply the givens rotations to a vector
%y and optionally a lower triangular matrix L.
[m,n] = size(R);
Rj = R(:,j);
R(:,j) = R(:,n);
R(:,n) = Rj;
if(nargin == 4)
    Lj = L(:,j);
    L(:,j) = L(:,n);
    L(:,n) = Lj;
end

for i=n-1:-1:j
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

for i = j+1:n-1
    r=sqrt(R(i,i)*R(i,i)+R(i+1,i)*R(i+1,i));
    c=R(i,i)/r;
    s=R(i+1,i)/r;
    G=[c,s; -s,c];
    R(i,i)=r;
    R(i+1,i)=0;
    R([i,i+1],i+1:n)=G*R([i,i+1],i+1:n);
    %addflops(flops_mul(G,R([i,i+1],i+1:n)));
    
    y([i,i+1]) = G*y([i,i+1]);
    %addflops(flops_mul(G,y([i,i+1])));
    
    if(nargin==4)
        L([i,i+1],[1:i,n]) = G*L([i,i+1],[1:i,n]);
        %addflops(flops_mul(G,L([i,i+1],[1:i,n])));
    end
end
