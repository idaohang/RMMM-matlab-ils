function [estimate, numExpanded,lastNode] = purebfsearchconstrained(R,y,lower,upper,maxNodes)
%
% [estimate,numExpanded,lastNode] = bfsearch(R,y,mem) produces the optimal solution to the upper triangular
%        integer least squares problem min_{z}||y-Rz|| by a search algorithm.
%
% Input arguments:
%    R ---- n by n real nonsingular upper triangular matrix
%    y ---- n-dimensional real vector
%    maxNodes ---- a parameter specifying the maximum number of nodes to
%                  keep in memory at any given time... must be set
%                  appropriately for given problem size
%
% Output arguments:
%    estimate - n by 1 integer vector (in double precision) specifying the optimal solution
%    numExpanded - An integer specifying the number of nodes visited during the
%    search
%    lastNode - An integer telling how much memory was required to complete
%    the search

    n = size(y,1);
    estimate = zeros(1,n);
    rad = inf; %set initial radius to infinity
    numExpanded = 1; %number of nodes expanded so far, 1 because we assume root expanded
    lastNode = 1; %last node in the node list

    PQb = pq_create(maxNodes); %Priority queue for currently active nodes, cost is partial distance to best child
    
    %Sigma is the lower bound from Chang's paper... it is the same for all
    %nodes at every level, so we can safely use it in the BFS!
    %Compute sigma for all values of k = 1:n - this is used in new upper bound
    %derived in the paper  
    sigmas = zeros(n,1);
    for i=2:n
        tempSum1 = 0;
        tempSum2 = 0;
        for j = i-1:n
            tempSum1 = tempSum1 + min(R(i-1,j)*lower,R(i-1,j)*upper);
            tempSum2 = tempSum2 + max(R(i-1,j)*lower,R(i-1,j)*upper);
        end
        if(sign((y(i-1) - tempSum1)) == sign((y(i-1) - tempSum2)))
            sigmas(i) = sigmas(i-1) + min((y(i-1) - tempSum1)^2,(y(i-1) - tempSum2)^2);
        end
    end
    
    intervalMid = zeros(1,maxNodes); %the middle of the valid integer range for a node (from SE algorithm)
    currentChild = zeros(1,maxNodes); %the current increment away from the mid that this node is at. It will follow the sequence 0,+1,-1,+2,-2... or 0,-1,+1,-2,+2...
    plusOrMinus = zeros(1,maxNodes); %records whether the first child was +1 (indicated by a positive number) or a -1 (indicated by a negative number)
    treeLevel = zeros(1,maxNodes); %the level in the tree that this node sits at
    cost = zeros(1,maxNodes); %the accumulated cost of this node. It is equal to the cost of the parent plus the cost to generate it
    nextChildCost = zeros(1,maxNodes); %the cost to generate the next child node
    symbolFromParent = zeros(1,maxNodes); %the symbol that was taken from the parent node to get to the current node
    parent = zeros(1,maxNodes); %the index of the nodes parent (0 for root)
    ubound = zeros(1,maxNodes); %vector to record whether we have hit the upper for this node
    lbound = zeros(1,maxNodes); %vector to record whether we hit lower for this node
    
    %Initialize the root node
    temp = y(n)/R(n,n);
    intervalMid(1) = round(temp);
    if(intervalMid(1) >= upper)
        ubound(1) = 1;
        intervalMid(1) = upper;
        plusOrMinus(1) = -1;
    else
        if(intervalMid(1) <= lower)
            lbound(1) = 1;
            intervalMid(1) = lower;
            plusOrMinus(1) = 1;
        else
            plusOrMinus(1) = temp - round(temp);
        end
    end
    currentChild(1) = 0;
    treeLevel(1) = n+1; %tree levels are one higher than matrix indices
    cost(1) = 0;
    nextChildCost(1) = (R(n,n)*intervalMid(1)-y(n))^2;
    symbolFromParent(1) = 0;
    parent(1) = 0;
    pq_push(PQb,1,nextChildCost(1));

    while(true)      
        %Pop from B into S until the node that comes out has a level less
        %than or equal to the maximum
        curNode = pq_pop(PQb);
        if(numExpanded == maxNodes)
            intervalMid(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %the middle of the valid integer range for a node (from SE algorithm)
            currentChild(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %the current increment away from the mid that this node is at. It will follow the sequence 0,+1,-1,+2,-2... or 0,-1,+1,-2,+2...
            plusOrMinus(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %records whether the first child was +1 (indicated by a positive number) or a -1 (indicated by a negative number)
            treeLevel(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %the level in the tree that this node sits at
            cost(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %the accumulated cost of this node. It is equal to the cost of the parent plus the cost to generate it
            nextChildCost(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %the cost to generate the next child node
            symbolFromParent(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %the symbol that was taken from the parent node to get to the current node
            parent(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %the index of the nodes parent (0 for root)
            ubound(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %vector to record whether we have hit the upper for this node
            lbound(maxNodes+1:2*maxNodes) = zeros(1,maxNodes); %vector to record whether we hit lower for this node
            maxNodes = 2*maxNodes;
        end
        if(nextChildCost(curNode) < rad) %if the cost to generate the next child stays within our radius
           if(treeLevel(curNode) > 2) %if the next best child is not a leaf, calculate its attributes and push it to PQb
                
                numExpanded = numExpanded+1;
                curChildIdx = numExpanded;
                
                
                currentChild(curChildIdx) = 0;
                treeLevel(curChildIdx) = treeLevel(curNode)-1;
                i = treeLevel(curChildIdx) - 1; %tree levels are one higher than matrix indices...
                cost(curChildIdx) = nextChildCost(curNode);
                symbolFromParent(curChildIdx) = currentChild(curNode) + intervalMid(curNode);                
                parent(curChildIdx) = curNode;
                
                if(ubound(curNode) ~= 1 || lbound(curNode) ~= 1)
                    if(lbound(curNode) == 1)
                        if(currentChild(curNode) + intervalMid(curNode) == upper)
                            ubound(curNode) = 1;
                        else
                            currentChild(curNode) = currentChild(curNode) + 1;
                        end
                    else
                        if(ubound(curNode) == 1)
                            if(currentChild(curNode) + intervalMid(curNode) == lower)
                                lbound(curNode) = 1;
                            else
                                currentChild(curNode) = currentChild(curNode) - 1; 
                            end
                        else
                            if(currentChild(curNode) == 0)
                                if(plusOrMinus(curNode) <= 0)
                                    currentChild(curNode) = -1;
                                    if(currentChild(curNode) + intervalMid(curNode) < lower)
                                        lbound(curNode) = 1;
                                        currentChild(curNode) = 1;
                                    end
                                else
                                    currentChild(curNode) = 1;
                                    if(currentChild(curNode) + intervalMid(curNode) > upper)
                                        ubound(curNode) = 1;
                                        currentChild(curNode) = -1;
                                    end
                                end
                            else
                                if(currentChild(curNode) < 0)
                                    if(plusOrMinus(curNode) <= 0)
                                        currentChild(curNode) = -currentChild(curNode);
                                        if(currentChild(curNode) + intervalMid(curNode) > upper )
                                            ubound(curNode) = 1;
                                            currentChild(curNode) = -currentChild(curNode) - 1;
                                        end
                                    else
                                        currentChild(curNode) = -currentChild(curNode) + 1;    
                                        if(currentChild(curNode) + intervalMid(curNode) > upper )
                                            ubound(curNode) = 1;
                                            currentChild(curNode) = -currentChild(curNode);
                                        end
                                    end
                                else
                                    if(plusOrMinus(curNode) <= 0)
                                        currentChild(curNode) = -currentChild(curNode) - 1;
                                        if(currentChild(curNode) + intervalMid(curNode) < lower )
                                            lbound(curNode) = 1;
                                            currentChild(curNode) = -currentChild(curNode);
                                        end
                                    else
                                        currentChild(curNode) = -currentChild(curNode);
                                        if(currentChild(curNode) + intervalMid(curNode) < lower )
                                            lbound(curNode) = 1;
                                            currentChild(curNode) = -currentChild(curNode) + 1;
                                        end
                                    end 
                                end
                            end 
                        end
                    end

                    if(lbound(curNode) ~= 1 || ubound(curNode) ~= 1)
                        %Calculate the cost of the new next child for the current
                        %node 
                        i=i+1; %needs to be increased temporarily because we are working at the parent level - 'i' was set to the childs level
                        sum=R(i,i)*(intervalMid(curNode) + currentChild(curNode));
                        curSymbolNode = curNode;
                        for j=i+1:n
                            sum = sum + R(i,j)*symbolFromParent(curSymbolNode);
                            curSymbolNode = parent(curSymbolNode);
                        end
                        temp = cost(curNode) + (sum - y(i))^2;
                        nextChildCost(curNode) = temp;


                        pq_push(PQb,curNode,nextChildCost(curNode)+sigmas(i)); %we push curnode back into the PQ with its new cost
                        i=i-1;
                    end
                end
                %Calculate the midpoint of the interval for the new child
                sum = 0;
                curSymbolNode = curChildIdx;
                for j=i+1:n
                    sum = sum+R(i,j)*symbolFromParent(curSymbolNode);
                    curSymbolNode = parent(curSymbolNode);
                end
                temp = (y(i) - sum)/R(i,i);
                intervalMid(curChildIdx) = round(temp);
                if(intervalMid(curChildIdx) >= upper)
                    ubound(curChildIdx) = 1;
                    intervalMid(curChildIdx) = upper;
                    plusOrMinus(curChildIdx) = -1;
                else
                    if(intervalMid(curChildIdx) <= lower)
                        lbound(curChildIdx) = 1;
                        intervalMid(curChildIdx) = lower;
                        plusOrMinus(curChildIdx) = 1;
                    else
                        plusOrMinus(curChildIdx) = temp - round(temp);
                    end
                end
                
                %Calculate the next child cost for the new nodes best
                %child - the sum above is re-used but with an extra term,
                %so just add it in
                sum = sum + R(i,i)*intervalMid(curChildIdx);                
                nextChildCost(curChildIdx) = cost(curChildIdx) + (sum - y(i))^2;
                pq_push(PQb,curChildIdx,nextChildCost(curChildIdx)+sigmas(i));
                                         
           else %we have hit a leaf... backtrack to find the resulting vector
                estimate(1) = intervalMid(curNode) + currentChild(curNode);
                %trace up the tree to get the current estimate
                for i=2:n
                    estimate(i) = symbolFromParent(curNode);
                    curNode = parent(curNode);
                end
                pq_delete(PQb);
                return;
           end
        else %cost to generate the next best child in PQb does not stay under current sphere radius... trim the current node
            continue;
        end
        
    end
    pq_delete(PQb);
end