function ephd=yuma2mat()
% DESCRIPTION: convert raw ephemeris data
%              into an internal ephemeris matrix, ephd(i,j).
% AUTHOR:      Valeri Perepetchai 12/08/99.

format short e
fid=fopen('austria.gpy','r');
if fid==-1,
 error('Error in opening the input yuma file')
end;
j = 1;
m=[];
while (feof(fid) == 0)
 s=fgets(fid);
 if ((s(1) ~= ' ') & (s(1) ~= '\n') & (s(1) ~= '*'))
  for i = 1:length(s)
   if (s(i) == ':')
    ss = s((i+1):length(s));
    m(j) = str2num(ss);
    j = j + 1;
    break
   end
  end
 end
end
fclose(fid);
n=(j-1)/13; %How many SVs?
ephd=zeros(n,13);
for k=1:n
ephd(k,1:13)=m((13*(k-1)+1):(k*13));
end

% Delete unhealthy data
i=1;
for k=1:n
    if(ephd(i,2)~=0)
       ephd=[ephd(1:(i-1),:);ephd((i+1):n,:)];
       n=n-1;
       i=i-1;
    end
    i=i+1; 
end
