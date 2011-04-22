runs = 1000;
time1 = zeros(1,runs);
time2 = zeros(1,runs);
expand1 = zeros(1,runs);
expand2 = zeros(1,runs);
sz = 40;
epochs = 100;
for i = 1:runs
    i
    [time1(i),expand1(i)] = frame(sz,epochs,0);
    [time2(i),expand2(i)] = frame(sz,epochs,1);
end