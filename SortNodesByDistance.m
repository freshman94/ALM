function [sortedNodes] = SortNodesByDistance(NotIn,targetNodePos,ClusterMatrix_,...
                        ClusterIdx,nodes)
%%功能：将nodes中的节点按照与targeNode的距离升序排序
%输入：
%NotIn用于指示nodes是否包含targeNodePos，1表示不包含
%nodes表示要加入的簇原本的结点序号
%输出：
%sortedNodes:2维矩阵，包含节点的序号和与结点targeNode的距离
%     fprintf('======================SortNodesByDistance====================\n');
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    nodeNum = size(nodes,2);
    %不包括targeNode的结点数量
    otherNodesNum = nodeNum+NotIn-1;
    nodesWithD = [];
    for i = 1:nodeNum
        if NotIn == 1 || targetNodePos(1) ~= nodes(i)
            distance = ( (targetNodePos(2)-Cluster(2,nodes(i)))^2 +...
                        (targetNodePos(3)-Cluster(3,nodes(i)))^2)^0.5;
            nodesWithD = [nodesWithD,[nodes(i),distance]'];
        end
    end
    
    sortedNodes = [];
    if otherNodesNum > 0
        distanceRow = nodesWithD(2,:);
        sortedNodesWithD = sort(distanceRow);
        sortedNodes = zeros(2,otherNodesNum);
        for i = 1:otherNodesNum
            pos = find(distanceRow==sortedNodesWithD(i));
            sortedNodes(:,i) = nodesWithD(:,pos);
        end
    end
end