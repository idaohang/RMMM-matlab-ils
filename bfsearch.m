function [estimate, numExpanded,lastNode] = bfsearch(R,y,mem,maxNodes,test,test2)
%
% [estimate,numExpanded,lastNode] = bfsearch(R,y,mem) produces the optimal solution to the upper triangular
%        integer least squares problem min_{z}||y-Rz|| by a search algorithm.
%
% Input arguments:
%    R ---- n by n real nonsingular upper triangular matrix
%    y ---- n-dimensional real vector
%    mem ---- a parameter which specifies the maximum amount of memory to
%             use for storing partial candidate solutions. Higher gives better
%             performance but requires more memory.
%    maxNodes ---- a parameter specifying the maximum number of nodes to
%                  keep in memory at any given time... must be set appropriately for a
%                  given mem or the program will crash
%    test ---- if true, when nodes are transferred from S to B, we will
%           transfer all nodes in the uppermost level. If false, we only transfer
%           the first 'mem' nodes from S to B where if there is a tie at the last
%           level, all the nodes after the 'mem'th are just left in S
%    test2 ---- if true then parents are only inserted back into PQb if
%               they still remain in the bottommost 'mem' nodes. if false
%               they are always put back in B. The advantage of always
%               putting them back is we may prune more nodes later.
%
% Output arguments:
%    estimate - n by 1 integer vector (in double precision) specifying the optimal solution
%    numExpanded - An integer specifying the number of nodes visited during the
%    search
%    lastNode - An integer telling how much memory was required to complete
%    the search

    n = size(y,1);
    estimate = zeros(1,n);
    maxActiveSize = mem; %this parameter controls the memory usage by limiting the next node to be expanded to lower levels. This forces the search to move down the tree. When a leaf is hit we can prune and tighten radius.
    rad = inf; %set initial radius to infinity
    numExpanded = 1; %number of nodes expanded so far, 1 because we assume root expanded
    lastNode = 1; %last node in the node list
    
    S = zeros(1,maxNodes); %List of currently inactive nodes (nodes not being considered for expansion currently)
    sPosition = 0; %current position in list S
    PQb = pq_create(maxNodes); %Priority queue for currently active nodes, cost is partial distance to best child
    nodesPerLevel = zeros(1,n+1); %array to keep track of the nodes at each level in the PQb
    
    freeMemory = zeros(1,maxNodes); %Array to keep track of free memory positions in the node array
    freeSize = 0; %current number of elements in the free memory array, we will add and remove from the back
    
    intervalMid = zeros(1,maxNodes); %the middle of the valid integer range for a node (from SE algorithm)
    currentChild = zeros(1,maxNodes); %the current increment away from the mid that this node is at. It will follow the sequence 0,+1,-1,+2,-2... or 0,-1,+1,-2,+2...
    plusOrMinus = zeros(1,maxNodes); %records whether the first child was +1 (indicated by a positive number) or a -1 (indicated by a negative number)
    treeLevel = zeros(1,maxNodes); %the level in the tree that this node sits at
    cost = zeros(1,maxNodes); %the accumulated cost of this node. It is equal to the cost of the parent plus the cost to generate it
    nextChildCost = zeros(1,maxNodes); %the cost to generate the next child node
    symbolFromParent = zeros(1,maxNodes); %the symbol that was taken from the parent node to get to the current node
    parent = zeros(1,maxNodes); %the index of the nodes parent (0 for root)
    inB = zeros(1,maxNodes);
    
    %Initialize the root node
    temp = y(n)/R(n,n);
    plusOrMinus(1) = temp - round(temp);
    intervalMid(1) = round(temp);
    currentChild(1) = 0;
    treeLevel(1) = n+1; %tree levels are one higher than matrix indices
    cost(1) = 0;
    nextChildCost(1) = (R(n,n)*intervalMid(1)-y(n))^2;
    symbolFromParent(1) = 0;
    parent(1) = 0;
    pq_push(PQb,1,nextChildCost(1));
    nodesPerLevel(n+1) = 1;
    inB(1) = 1;
    
    while(pq_size(PQb) ~= 0)
        %Find the lowest k nodes in the PQb... find which level is the
        %maximum. The value in k after this loop represents the maximum
        %level from which we may select a node to be expanded.
        nodeCount=0;
        for k=1:n+1
            nodeCount = nodeCount + nodesPerLevel(k);
            if(nodeCount >= maxActiveSize)
                break;
            end
        end
        
        %Pop from B into S until the node that comes out has a level less
        %than or equal to the maximum
        curNode = pq_pop(PQb);
        inB(curNode) = 0;
        nodesPerLevel(treeLevel(curNode)) = nodesPerLevel(treeLevel(curNode))-1;
        
        while(treeLevel(curNode) > k)
            %before we add a node to the list S, we must first make sure
            %none of its ancestors are still in B. If they are, and B is
            %destroyed, this node will be corrupted.
            curSymbolNode = curNode;
            while(parent(curSymbolNode) ~= 0 && inB(parent(curSymbolNode)) == 1)
                pq_push(PQb,parent(curSymbolNode),0); %update the parents cost in the PQ to 0 so we can pop it out
                pq_pop(PQb); %remove the parent
                inB(parent(curSymbolNode)) = 0;
                sPosition = sPosition+1;
                S(sPosition) = parent(curSymbolNode); %add the parent node to the list S of currently un-used nodes
                nodesPerLevel(treeLevel(parent(curSymbolNode))) = nodesPerLevel(treeLevel(parent(curSymbolNode)))-1;
                curSymbolNode = parent(curSymbolNode); %move up in the tree
            end
            
            sPosition = sPosition+1;
            S(sPosition) = curNode; %add the current node to the list S of currently un-used nodes
            
            curNode = pq_pop(PQb); %retrieve the next candidate node for expansion
            inB(curNode) = 0;
            nodesPerLevel(treeLevel(curNode)) = nodesPerLevel(treeLevel(curNode))-1;
        end
        
        if(nextChildCost(curNode) < rad) %if the cost to generate the next child stays within our radius
           if(treeLevel(curNode) > 2) %if the next best child is not a leaf, calculate its attributes and push it to PQb
                numExpanded = numExpanded+1;
                
                %if there are no free memory locations expand the node
                %array
                if(freeSize == 0)
                    curChildIdx = lastNode + 1;
                    lastNode = lastNode+1;
                else %if there are free memory locations, use them
                    curChildIdx = freeMemory(freeSize);
                    freeSize = freeSize-1;
                end
                
                currentChild(curChildIdx) = 0;
                treeLevel(curChildIdx) = treeLevel(curNode)-1;
                i = treeLevel(curChildIdx) - 1; %tree levels are one higher than matrix indices...
                cost(curChildIdx) = nextChildCost(curNode);
                symbolFromParent(curChildIdx) = currentChild(curNode) + intervalMid(curNode);
                parent(curChildIdx) = curNode;
                
                %Update the parent nodes next best child and update the PQ
                %use the plusOrMinus value to determine which child should
                %be generated next
                if(currentChild(curNode) == 0)
                    if(plusOrMinus(curNode) <= 0)
                        currentChild(curNode) = -1;
                    else
                        currentChild(curNode) = 1;
                    end
                else
                    if(currentChild(curNode) < 0)
                        if(plusOrMinus(curNode) <= 0)
                            currentChild(curNode) = -currentChild(curNode);
                        else
                            currentChild(curNode) = -currentChild(curNode) + 1;                   
                        end
                    else
                        if(plusOrMinus(curNode) <= 0)
                            currentChild(curNode) = -currentChild(curNode) - 1;
                        else
                            currentChild(curNode) = -currentChild(curNode);                   
                        end 
                    end
                end
                
                %Calculate the cost of the new next child for the current node  
                %this calculation could be done MUCH more efficiently (no
                %loop) since only one term in the sum changes for the new
                %child. This could be a future improvement but the real
                %bottleneck is the priority queue.
                i=i+1; %needs to be increased temporarily because we are working at the parent level - 'i' was set to the childs level
                sum=R(i,i)*(intervalMid(curNode) + currentChild(curNode));
                curSymbolNode = curNode;
                for j=i+1:n
                    sum = sum + R(i,j)*symbolFromParent(curSymbolNode);
                    curSymbolNode = parent(curSymbolNode);
                end
                temp = cost(curNode) + (sum - y(i))^2;
                nextChildCost(curNode) = temp;
                
                if(test2)
                    if((nodeCount - nodesPerLevel(k) + 1 < maxActiveSize)) %if curnode is still in bottommost nodes, add it back to B
                        pq_push(PQb,curNode,nextChildCost(curNode)); %we push curnode back into the PQ... we could have also just never popped it and changed the cost with an update but the cost would be similar
                        inB(curNode) = 1;
                        nodesPerLevel(treeLevel(curNode)) = nodesPerLevel(treeLevel(curNode))+1;
                    else %put curNode into S
                        %before we add a node to the list S, we must first make sure
                        %none of its ancestors are still in B. If they are, and B is
                        %destroyed, this node will be corrupted.
                        curSymbolNode = curNode;
                        while(parent(curSymbolNode) ~= 0 && inB(parent(curSymbolNode)) == 1)
                            pq_push(PQb,parent(curSymbolNode),0); %update the parents cost in the PQ to 0 so we can pop it out
                            pq_pop(PQb); %remove the parent
                            inB(parent(curSymbolNode)) = 0;
                            sPosition = sPosition+1;
                            S(sPosition) = parent(curSymbolNode); %add the parent node to the list S of currently un-used nodes
                            nodesPerLevel(treeLevel(parent(curSymbolNode))) = nodesPerLevel(treeLevel(parent(curSymbolNode)))-1;
                            curSymbolNode = parent(curSymbolNode); %move up in the tree
                        end

                        sPosition = sPosition+1;
                        S(sPosition) = curNode; %add the current node to the list S of currently un-used nodes
                    end
                else
                    pq_push(PQb,curNode,nextChildCost(curNode)); %we push curnode back into the PQ... we could have also just never popped it and changed the cost with an update but the cost would be similar
                    inB(curNode) = 1;
                    nodesPerLevel(treeLevel(curNode)) = nodesPerLevel(treeLevel(curNode))+1;
                end
                
                i=i-1;
                
                %Calculate the midpoint of the interval for the new child
                sum = 0;
                curSymbolNode = curChildIdx;
                for j=i+1:n
                    sum = sum+R(i,j)*symbolFromParent(curSymbolNode);
                    curSymbolNode = parent(curSymbolNode);
                end
                temp = (y(i) - sum)/R(i,i);
                intervalMid(curChildIdx) = round(temp);
                plusOrMinus(curChildIdx) = temp - intervalMid(curChildIdx);
                
                %Calculate the next child cost for the new nodes best
                %child - the sum above is re-used but with an extra term,
                %so just add it in
                sum = sum + R(i,i)*intervalMid(curChildIdx);                
                nextChildCost(curChildIdx) = cost(curChildIdx) + (sum - y(i))^2;
                nodesPerLevel(i+1) = nodesPerLevel(i+1)+1;
                pq_push(PQb,curChildIdx,nextChildCost(curChildIdx));
                inB(curChildIdx) = 1;
                                         
           else %we have hit a leaf... backtrack to find the resulting vector, update the radius and refresh B
                estimate(1) = intervalMid(curNode) + currentChild(curNode);
                tempRad = nextChildCost(curNode);
                if (tempRad < rad)
                    rad = tempRad;
                end
                %trace up the tree to get the current estimate
                for i=2:n
                    estimate(i) = symbolFromParent(curNode);
                    curNode = parent(curNode);
                end
                %remove all nodes in B from search
                %first set all indices currently in B to be available as
                %free memory
                freeMemory(freeSize+1:freeSize+pq_size(PQb)) = pq_getArray(PQb);
                freeSize = freeSize+pq_size(PQb);
                pq_delete(PQb);
                nodesPerLevel = zeros(1,n+1);
                inB = zeros(1,maxNodes);
                PQb = pq_create(maxNodes);
                
                %Sort S by tree level so if we take a chunk of nodes of size
                %maxActiveSize from the end of S, we get the bottom-most
                %'maxActiveSize' nodes in the tree
                [~,ind] = sort(treeLevel(S(1:sPosition)),'descend');
                S = S(ind);
                start = max(1,sPosition-maxActiveSize);

                %It is arbitrary to pick the last nodes because there are many
                %ties in the tree level. So, for any other nodes in S that have
                %an appropriate tree level, include them too
                if(test)
                    while(start > 1 && treeLevel(S(start)) == treeLevel(S(start-1)))
                        start = start - 1;
                    end
                end
                
                %Now empty the appropriate nodes from S to B
                for i=start:sPosition
                    curNode = S(i);
                    pq_push(PQb,curNode,nextChildCost(curNode));
                    inB(curNode) = 1;
                    nodesPerLevel(treeLevel(curNode)) = nodesPerLevel(treeLevel(curNode))+1;
                end
                sPosition = start - 1;
           end
        else %cost to generate the next best child in PQb does not stay under current sphere radius... throw away and re-populate B 
            %first set all indices currently in B to be available as
            %free memory
            freeMemory(freeSize+1:freeSize+pq_size(PQb)) = pq_getArray(PQb);
            freeSize = freeSize+pq_size(PQb); 
            pq_delete(PQb);
            nodesPerLevel = zeros(1,n+1);
            inB = zeros(1,maxNodes);
            PQb = pq_create(maxNodes);
           
            %Sort S by tree level so if we take a chunk of nodes of size
            %maxActiveSize from the end of S, we get the bottom-most
            %'maxActiveSize' nodes in the tree
            [~,ind] = sort(treeLevel(S(1:sPosition)),'descend');
            S = S(ind);
            start = max(1,sPosition-maxActiveSize);
       
            %It is arbitrary to pick the last nodes because there are many
            %ties in the tree level. So, for any other nodes in S that have
            %an appropriate tree level, include them too
            if(test)
                while(start > 1 && treeLevel(S(start)) == treeLevel(S(start-1)))
                    start = start - 1;
                end
            end
            
            %Now move the appropriate ndoes from S to B
            for i=start:sPosition
                curNode = S(i);
                pq_push(PQb,curNode,nextChildCost(curNode));
                inB(curNode) = 1;
                nodesPerLevel(treeLevel(curNode)) = nodesPerLevel(treeLevel(curNode))+1;
            end
            sPosition = start - 1;

        end
        
    end
    
end