function [c,s,r]=given(a,b)
if a==0
   c=1;s=0;r=b;
else
   r=sqrt(a*a+b*b);
   s=-a/r;
   c=b/r;
end
