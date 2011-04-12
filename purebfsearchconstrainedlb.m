function [estimate, numExpanded,lastNode] = purebfsearchconstrainedlb(R,y,lower,upper,maxNodes)
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
    m = n;
    estimate = zeros(1,n);
    rad = inf; %set initial radius to infinity
    numExpanded = 1; %number of nodes expanded so far, 1 because we assume root expanded
    lastNode = 1; %last node in the node list

    PQb = pq_create(maxNodes); %Priority queue for currently active nodes, cost is partial distance to best child
    
    %Preprocessing for lower bound
    [eVals F] = preProcessLB(R);
    tempF = zeros(1,m);
    
    nextChildF = zeros(maxNodes,m); %a value which we need to propagate to efficiently compute the lower-bounds - this is for the next child
    nodeF = zeros(maxNodes,m); %same as above but for the node itself, not its child
    lowerb = zeros(1,maxNodes); %store the LB for this node to use for its children
    nextChildLB = zeros(1,maxNodes); %store the lb for the next child
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
    parent(1) = 1;
    
    %Calculate the lowerbound for the roots first child
    nextChildF(1,1:m-1) = F(1:m-1,1:m-1)*(y(1:m-1) - R(1:m-1,m)*intervalMid(1));
    nodeF(1,1:m) = zeros(1,m);
    for i = 1:m-1
        tempF(i) = min(round(nextChildF(1,i)),upper);
        tempF(i) = max(tempF(i),lower); 
    end
    nextChildLB(1) = eVals(m-1)*norm((tempF(1:m-1) - nextChildF(1,1:m-1)))^2;
    
    pq_push(PQb,1,nextChildCost(1));
    
    while(true)
       
        %Pop from B into S until the node that comes out has a level less
        %than or equal to the maximum
        curNode = pq_pop(PQb);
        
        if(nextChildCost(curNode) < rad) %if the cost to generate the next child stays within our radius
           if(treeLevel(curNode) > 2) %if the next best child is not a leaf, calculate its attributes and push it to PQb
                
                numExpanded = numExpanded+1;
                curChildIdx = numExpanded;
                
                currentChild(curChildIdx) = 0;
                treeLevel(curChildIdx) = treeLevel(curNode)-1;
                i = treeLevel(curChildIdx) - 1; %tree levels are one higher than matrix indices...
                cost(curChildIdx) = nextChildCost(curNode);
                nodeF(curChildIdx,1:i) = nextChildF(curNode,1:i);
                symbolFromParent(curChildIdx) = currentChild(curNode) + intervalMid(curNode);
                parent(curChildIdx) = curNode;
                lowerb(curChildIdx) = nextChildLB(curNode);
                
                if(ubound(curNode) ~= 1 || lbound(curNode) ~= 1)
                    if(lbound(curNode) == 1)
                        currentChild(curNode) = currentChild(curNode) + 1;
                        if(currentChild(curNode) + intervalMid(curNode) == upper)
                            ubound(curNode) = 1;
                        end
                    else
                        if(ubound(curNode) == 1)
                            currentChild(curNode) = currentChild(curNode) - 1; 
                            if(currentChild(curNode) + intervalMid(curNode) == lower)
                                lbound(curNode) = 1;
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
                                            currentChild(curNode) = -currentChild(curNode) - 1;
                                        end
                                    end 
                                end
                            end 
                        end
                    end

                    
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
                    
                    %Calculate the lowerbound for the next child                           
                    lb = 0;
                    if(i == m || i == 1)
                        if(i == m)
                            nextChildF(curNode,1:i-1) = F(1:m-1,1:m-1)*(y(1:m-1) - R(1:m-1,m)*(currentChild(curNode) + intervalMid(curNode)));
                            for p = 1:i-1
                                tempF(p) = min(round(nextChildF(curNode,p)),upper);
                                tempF(p) = max(tempF(p),lower); 
                            end
                            nextChildLB(curNode) = eVals(i-1)*norm((tempF(1:i-1) - nextChildF(curNode,1:i-1)))^2;
                        end
                    else
                        nextChildF(curNode,1:i-1) = nodeF(curNode,1:i-1)' - F(1:i-1,i)*(y(i) - sum );
                        for p = 1:i-1
                            tempF(p) = min(round(nextChildF(curNode,p)),upper);
                            tempF(p) = max(tempF(p),lower); 
                        end
                        nextChildLB(curNode) = eVals(i-1)*norm((tempF(1:i-1) - nextChildF(curNode,1:i-1)))^2;
                    end

                    pq_push(PQb,curNode,nextChildCost(curNode)+lowerb(parent(curNode))); %we push curnode back into the PQ with its new cost
                    i=i-1;
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

                %Calculate the LB for the new child's next child
                sum = sum + R(i,i)*intervalMid(curChildIdx);                

                lb = 0;
                if(i == m || i == 1)
                    if(i == m)
                        nextChildF(curChildIdx,1:i-1) = F(1:m-1,1:m-1)*(y(1:m-1) - R(1:m-1,m)*(intervalMid(curChildIdx)));
                        for p = 1:i-1
                            tempF(p) = min(round(nextChildF(curChildIdx,p)),upper);
                            tempF(p) = max(tempF(p),lower); 
                        end
                        nextChildLB(curChildIdx) = eVals(i-1)*norm((tempF(1:i-1) - nextChildF(curChildIdx,1:i-1)))^2;
                    end
                else
                    nextChildF(curChildIdx,1:i-1) = nodeF(curChildIdx,1:i-1)' - F(1:i-1,i)*(y(i) -sum);
                    for p = 1:i-1
                        tempF(p) = min(round(nextChildF(curChildIdx,p)),upper);
                        tempF(p) = max(tempF(p),lower); 
                    end
                    nextChildLB(curChildIdx) = eVals(i-1)*norm((tempF(1:i-1) - nextChildF(curChildIdx,1:i-1)))^2;
                end
                
                %Calculate the next child cost for the new nodes best
                %child - the sum above is re-used but with an extra term,
                %so just add it in
                nextChildCost(curChildIdx) = cost(curChildIdx) + (sum - y(i))^2;
                pq_push(PQb,curChildIdx,nextChildCost(curChildIdx)+lowerb(parent(curChildIdx)));
                                         
           else %we have hit a leaf... backtrack to find the resulting vector
                estimate(1) = intervalMid(curNode) + currentChild(curNode);
                %trace up the tree to get the current estimate
                for i=2:n
                    estimate(i) = symbolFromParent(curNode);
                    curNode = parent(curNode);
                end
                return;
           end
        else %cost to generate the next best child in PQb does not stay under current sphere radius... trim the current node
            continue;
        end
        
    end
    
end