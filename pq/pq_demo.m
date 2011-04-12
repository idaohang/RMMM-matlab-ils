% PQ_DEMO compiles library and illustrate Priority Queue's functionalities
%
% Copyright (c) 2008 Andrea Tagliasacchi
% All Rights Reserved
% email: andrea.tagliasacchi@gmail.com
% $Revision: 1.1.1.1 $  Created on: May 22, 2009
mex pq/pq_create.cpp
mex pq/pq_delete.cpp
mex pq/pq_push.cpp
mex pq/pq_pop.cpp
mex pq/pq_size.cpp
mex pq/pq_top.cpp
mex pq/pq_getArray.cpp

clc, clear, close all;
pq = pq_create( 10 ); 
testVec = zeros(1,100);

for i=1:100
    disp(sprintf('\n')); %newline

    %--- create a random entry
    cost = rand(1);
    pq_push(pq, i, cost);
    testVec(i) = cost;
    fprintf('idx: %f cost: %f \n',i,cost);
 
end

[testVec idx] = sort(testVec);
for i=1:100
    if(idx(i) - pq_pop(pq) ~= 0)
        error('Wrong!');
    else
        fprintf('%f good\n',i);
    end
end
% disp(pq_size(pq));
% [idx cost] = pq_pop(pq);
% fprintf('idx: %f cost: %f \n',idx,cost);
% pq_push(pq,idx,cost);
% disp(pq_size(pq));

pq_delete(pq);