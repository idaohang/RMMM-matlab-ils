function [estimate, numExpanded,lastNode] = purebfsearchlb(R,y,maxNodes)
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
    
    intervalMid = zeros(1,maxNodes); %the middle of the valid integer range for a node (from SE algorithm)
    currentChild = zeros(1,maxNodes); %the current increment away from the mid that this node is at. It will follow the sequence 0,+1,-1,+2,-2... or 0,-1,+1,-2,+2...
    plusOrMinus = zeros(1,maxNodes); %records whether the first child was +1 (indicated by a positive number) or a -1 (indicated by a negative number)
    treeLevel = zeros(1,maxNodes); %the level in the tree that this node sits at
    cost = zeros(1,maxNodes); %the accumulated cost of this node. It is equal to the cost of the parent plus the cost to generate it
    nextChildCost = zeros(1,maxNodes); %the cost to generate the next child node
    nextChildF = zeros(maxNodes,m); %a value which we need to propagate to efficiently compute the lower-bounds - this is for the next child
    nodeF = zeros(maxNodes,m); %same as above but for the node itself, not its child
    symbolFromParent = zeros(1,maxNodes); %the symbol that was taken from the parent node to get to the current node
    parent = zeros(1,maxNodes); %the index of the nodes parent (0 for root)
    
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
    
    %Calculate the lowerbound for the roots first child
    nextChildF(1,1:m-1) = F(1:m-1,1:m-1)*(y(1:m-1) - R(1:m-1,m)*intervalMid(1));
    nodeF(1,1:m) = zeros(1,m);
    lb = eVals(m-1)*norm((round(nextChildF(1,1:m-1)) - nextChildF(1,1:m-1)))^2;
    
    pq_push(PQb,1,nextChildCost(1)+lb);
    
    while(true)
       
        %Pop from B into S until the node that comes out has a level less
        %than or equal to the maximum
        curNode = pq_pop(PQb);
        
        if(nextChildCost(curNode) < rad) %if the cost to generate the next child stays within our radius... pointless check for BFS! - radius is always infinite
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
                        lb = eVals(i-1)*norm((round(nextChildF(curNode,1:i-1)) - nextChildF(curNode,1:i-1)))^2;
                    end
                else
                    nextChildF(curNode,1:i-1) = nodeF(curNode,1:i-1)' - F(1:i-1,i)*(y(i) - sum );
                    lb = eVals(i-1)*norm((round(nextChildF(curNode,1:i-1)) - nextChildF(curNode,1:i-1)))^2;
                end
                
                %The lowerbound is added to the cost - A* algorithm
                pq_push(PQb,curNode,nextChildCost(curNode)+lb); %we push curnode back into the PQ with its new cost
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
         
                
                %Calculate the LB for the new child's next child  
                sum = sum + R(i,i)*intervalMid(curChildIdx);                

                lb = 0;
                if(i == m || i == 1)
                    if(i == m)
                        nextChildF(curChildIdx,1:i-1) = F(1:m-1,1:m-1)*(y(1:m-1) - R(1:m-1,m)*(intervalMid(curChildIdx)));
                        lb = eVals(i-1)*norm((round(nextChildF(curChildIdx,1:i-1)) - nextChildF(curChildIdx,1:i-1)))^2;
                    end
                else
                    nextChildF(curChildIdx,1:i-1) = nodeF(curChildIdx,1:i-1)' - F(1:i-1,i)*(y(i) - sum );
                    lb = eVals(i-1)*norm((round(nextChildF(curChildIdx,1:i-1)) - nextChildF(curChildIdx,1:i-1)))^2;
                end
               
                %Calculate the next child cost for the new nodes best
                %child - the sum above is re-used but with an extra term,
                %so just add it in
                nextChildCost(curChildIdx) = cost(curChildIdx) + (sum - y(i))^2;
                
                %Again push into the PQ the cost plus the lowerbound
                pq_push(PQb,curChildIdx,nextChildCost(curChildIdx)+lb);
                                         
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