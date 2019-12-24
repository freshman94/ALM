%%计算结点的优先值
function [priority] = GetPriority(nodeIdx,vertexStability,am,vertexDelay,interClusterInfo)

%     fprintf('====================GetPriority==============================\n');
%各个参数的权值
    a = 0.4;
    b = 0.4;
    c = 1 - a - b;
    
    stability = vertexStability(nodeIdx);
    stabilitySum = sum(vertexStability);
    nodesNum = size(vertexStability,2);
    
    %统计结点以及其所有相邻节点度之和，与所有节点的度之和进行归一化
    adjNodes = [];
    for i = 1:nodesNum
        if am(nodeIdx,i) == 1
            adjNodes = [adjNodes,i];
        end
    end
    adjNodesNum = size(adjNodes,2);
    degree = sum(am(nodeIdx,:));
    degree = degree + size(interClusterInfo{3,nodeIdx},2);
    for i = 1:adjNodesNum
        degree = degree + sum(am(adjNodes(i),:));
        degree = degree + size(interClusterInfo{3,adjNodes(i)},2);
    end
    degreeSum = 0;
    for i = 1:nodesNum
        degreeSum = degreeSum + sum(am(i,:));
        degreeSum = degreeSum + size(interClusterInfo{3,i},2);
    end
    
    delay = vertexDelay(nodeIdx);
    if isinf(delay)
        delay = 10e10;
    end
    delaySum = delay;
    for i = 1:nodesNum
%             fprintf('vertexDelay(%d) = %f\n',i,vertexDelay(i));
        delaySum = delaySum + vertexDelay(i);
    end
    
    priority = a*(1-stability/stabilitySum) + b*(1-degree/degreeSum) + ...
                c*delay/delaySum;
            
%     fprintf('nodeIdx = %d,stability = %f,stabilitySum = %f,degree = %d,degreeSum = %d,delay = %f,delaySum = %f,priority = %f\n',...
%         nodeIdx,stability,stabilitySum,degree,degreeSum,delay,delaySum,priority);
end