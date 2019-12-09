function [sortedNodes] = SortNodesByDistance(NotIn,targetNodePos,ClusterMatrix_,...
                        ClusterIdx,nodes)
%%���ܣ���nodes�еĽڵ㰴����targeNode�ľ�����������
%���룺
%NotIn����ָʾnodes�Ƿ����targeNodePos��1��ʾ������
%nodes��ʾҪ����Ĵ�ԭ���Ľ�����
%�����
%sortedNodes:2ά���󣬰����ڵ����ź�����targeNode�ľ���
%     fprintf('======================SortNodesByDistance====================\n');
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    nodeNum = size(nodes,2);
    %������targeNode�Ľ������
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