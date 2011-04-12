function [A,y,bx]=randILS(m,n,econd,mean,svar,snoise,x)

if nargin<2
    n=m;
end
if nargin<3
    econd=-1;
end
if nargin<4
    if econd <0
        mean=0;
    else
        mean=1;
    end
end
if nargin<5
    if econd<0
        svar=1;
    else
        svar=0.2;
    end
end
if nargin<6
    snoise=0.01;
end
if nargin<7
    x=round(20*rand(n,1))-10;
end

if econd<0
A=mean+svar*randn(m,n);
else
    mn=min(m,n);
    erange=linspace(-0.5*econd,0.5*econd,mn);
    brange=sort(1+9*rand(1,mn));
    %brange=1+9*rand(1,mn);
    rperm=randperm(mn);
    D=diag(brange(rperm).*(10.^erange(rperm)));
    %D=diag(brange.*(10.^erange));
    %maxD=abs(mean+svar*randn(1))*10^econd;
    %D=diag([1,maxD,rand(1,mn-2)*(maxD-1)+1]);
    [U,R]=qr(randn(m));
    [V,R]=qr(randn(n));
    A=U(:,1:mn)*D*V(1:mn,:);
end
y=A*x+snoise*randn(m,1);