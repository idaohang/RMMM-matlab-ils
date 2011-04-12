clear;
%Solve Test
k=4;
r = 16;
S = abs(eye(r) .* rand(r,r)*10) + 1;
q = rand(1,r)*10;

% Use the secant method to solve for a positive value of lamda
a = 0;
b = 0.5;
fa = 0;
fb = 0;
for i=1:r
    fa = fa + ((S(i,i)*q(i))/(S(i,i) + a))^2;
    fb = fb + ((S(i,i)*q(i))/(S(i,i) + b))^2;
end
fa = fa - k + 0.25;
fb = fb - k + 0.25;

count=0;
while(fb > 10^-8)
    count=count+1;
    if(a == b)
        break;
    end
    a = a - fa/((fb-fa)/(b-a));
    fa = 0;
    for i=1:r
        fa = fa + ((S(i,i)*q(i))/(S(i,i) + a))^2;
    end
    fa = fa - k + 0.25;

    if(a == b)
        break;
    end
    b = b - fb/((fa-fb)/(a-b));
    fb = 0;
    for i=1:r
        fb = fb + ((S(i,i)*q(i))/(S(i,i) + b))^2;
    end
    fb = fb - k + 0.25;           
end
count
fa=0;
for i=1:r
    fa = fa + ((S(i,i)*q(i))/(S(i,i) + b))^2;
end

fa - k + 0.25
