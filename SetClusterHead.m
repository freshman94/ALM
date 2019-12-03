function [IsClusterHead] = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority)
    clusterHead = IsClusterHead_{ClusterIdx(1),ClusterIdx(2)};
    %重置簇首
    nodesNum = size(vertexPriority,2);
    clusterHead = zeros(1,nodesNum);
    
    priority_sorted = sort(vertexPriority);
    %处理优先值为-inf的情况
    for i = 1:nodesNum
        pos = find(vertexPriority == priority_sorted(i));
        if isinf(vertexPriority(pos))
            continue;
        else
            clusterHead(pos) = 1;
            break;
        end
    end
    IsClusterHead_{ClusterIdx(1),ClusterIdx(2)} = clusterHead;
    IsClusterHead = IsClusterHead_;   
end