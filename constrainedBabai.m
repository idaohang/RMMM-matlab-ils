function [z] = constrainedBabai(R,y,lower,upper)
n = length(R);
% Initialization
z = zeros(n,1); % the current point
z(n) = max(min(round(y(n)/R(n,n)),upper),lower);

for i = n-1:-1:1
    z(i) = max(min(round((y(i) - R(i,i+1:n)*z(i+1:n))/R(i,i)),upper),lower);
end

end