function [z] = babai(R,y)
n = length(R);
% Initializes
z = zeros(n,1); % the current point
z(n) = round(y(n)/R(n,n));

for i = n-1:-1:1
    z(i) = round((y(i) - R(i,i+1:n)*z(i+1:n))/R(i,i));
end

end
