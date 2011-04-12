function [R,Z] = IGT(R,Z,j)

for i = j-1:-1:1
    alpha = round(R(i,j)/R(i,i));
    if alpha ~= 0
        R(1:i,j) = R(1:i,j) - alpha*R(1:i,i);
        Z(:,j) = Z(:,j) - alpha*Z(:,i);
    end
end
