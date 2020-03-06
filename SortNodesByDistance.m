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
    nodesWithD = zeros(otherNodesNum,2);
    j = 1;
    for i = 1:nodeNum
        if NotIn == 1 || targetNodePos(1) ~= nodes(i)
            distance = ( (targetNodePos(2)-Cluster(2,nodes(i)))^2 +...
                        (targetNodePos(3)-Cluster(3,nodes(i)))^2)^0.5;
            nodesWithD(j,:) = [nodes(i),distance];
            j = j + 1;
        end
    end
    
    sortedNodes = sortrows(nodesWithD,2);
    sortedNodes = sortedNodes';
    
%     sortedNodes = zeros(2,otherNodesNum);
%     if otherNodesNum > 0
%         distanceRow = nodesWithD(2,:);
%         sortedNodesWithD = sort(distanceRow);
%         sortedNodes = zeros(2,otherNodesNum);
%         for i = 1:otherNodesNum
%             pos = find(distanceRow==sortedNodesWithD(i));
%             if isempty(pos) || size(pos,2) ~= 1
%                 error('pos empty\n');
%             end
%             sortedNodes(:,i) = nodesWithD(:,pos);
%         end
%     end
end