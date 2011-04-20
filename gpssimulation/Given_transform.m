function R_bar = Given_transform(A)
global CC_m;
for i=4:3+CC_m-1
   for j=3:-1:1
      [c,s,r]=given(A(j,i),A(i,i));
      G=[c,s;-s,c];
      A([j,i],:) = G*A([j,i],:);
   end
end
R_bar = A(1:3,1:3);
