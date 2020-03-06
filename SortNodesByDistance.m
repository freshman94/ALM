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