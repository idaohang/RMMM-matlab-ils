function [Q,R,ind,y] = gradientReduction(B,y,estimate)
    grad = B'*(B*estimate - y); %Compute the gradient
    [~,ind] = sort(abs(grad),1,'ascend');
    B = B(:,ind);
    y = y(ind);
    [Q,R] = qr(B);
    y = Q'*y;
end