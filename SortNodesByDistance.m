function [sortedNodes] = SortNodesByDistance(NotIn,targetNodePos,ClusterMatrix_,ClusterIdx,nodes)
%% 将nodes中的节点按照与targeNode的距离升序排序
%NotIn用于指示targeNodePos是否有可能被nodes包含，1表示没有可能
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    nodeNum = size(nodes,2);
    nodesWithD = [];
    for i = 1:nodeNum
        if NotIn == 1 || targetNodePos(1) ~= nodes(i)
            distance = ( (targetNodePos(2)-Cluster(2,nodes(i)))^2 +...
                        (targetNodePos(3)-Cluster(3,nodes(i)))^2)^0.5;
            nodesWithD = [nodesWithD,[nodes(i),distance]'];
        end
    end
    sortedNodes = [];
    %nodes可能包含targeNode，因此需要重新size
    nodeNum = size(nodesWithD,2);
    if nodeNum > 0
        distanceRow = nodesWithD(2,:);
        sortedNodesWithD = sort(distanceRow);
        sortedNodes = zeros(1,nodeNum);
        for i = 1:nodeNum
            pos = find(distanceRow==sortedNodesWithD(i));
            sortedNodes(i) = nodesWithD(1,pos);
        end
    end
        
    %调试
    fprintf('============================SortNodesByDistance======================\n');
    for i = 1:nodeNum
        fprintf('%d\t',sortedNodes(i));
    end
    fprintf('\n===================================================================\n');
end