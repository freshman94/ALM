function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = VertexOutAndIn(MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexDelay_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%%输入参数
%ClusterIdx: 1*2行向量，指明某个簇
%nodeIdx: 数值，指明簇中某个节点

    ClusterOut = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    nodesOutNum = size(ClusterOut,2);
    %统计受影响节点
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    %判断受影响节点是否需要重新加入簇
    affectedNodes = [];
    for i = 1:nodesOutNum
        if am(nodeIdx,i) == 1 
            affectedNodes = [affectedNodes,i];
        end
    end  
    
    affectedNum = size(affectedNodes,2);
    
    %调试
    fprintf('clusterIdx=[%d,%d],nodeIdx = %d,affectedNum = %d\n',ClusterIdx(1),ClusterIdx(2),nodeIdx,affectedNum);
    
    fprintf('================print am==========\n');
    for i = 1:nodesOutNum
        for j = 1:nodesOutNum
            fprintf('am(%d,%d) = %d\t',i,j,am(i,j));
        end
        fprintf('\n');
    end
    fprintf('=========================\n');
    
    am(:,nodeIdx) = zeros(nodesOutNum,1);
    am(nodeIdx,:) = zeros(1,nodesOutNum);
    AM_{ClusterIdx(1),ClusterIdx(2)} = am;
    
    AllConnected = 1;
    for i = 1:affectedNum
        am = AM_{ClusterIdx(1),ClusterIdx(2)};
        %若节点的度为1，则其退出不会对其它结点有影响
        isConnected = CheckConnected(am,affectedNodes(i),nodeIdx);
        %调试代码
        fprintf('affectedNode = %d,isConnected=%d\n',affectedNodes(i),isConnected);
        %此变量用于指示脱离的节点是否重新加入成功
        %该节点脱离，重新加入
        if isConnected == 0
            affectedNodePos = ClusterOut(:,affectedNodes(i));
            [AllConnected,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_] = ...
            OldNodeJoin(MaxLinkDistance,affectedNodes,ClusterMatrix_,ClusterIdx,affectedNodes(i),affectedNodePos,nodeIdx,AM_,...
            EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB);
            %拓扑连通性重建成功
%             if AllConnected == 1
%                 break;
%             end
        end
    end
    
    fprintf('AllConnected = %d\n',AllConnected);
    
    %脱离的节点没有加入成功，则将nodeIdx的所有邻节点相连，以重新构建拓扑连通性
    if AllConnected == 0
        [AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_] = ...
        ConnectNodeIdxAllAdj(MaxLinkDistance,ClusterMatrix_,ClusterIdx,affectedNodes,AM_,...
            EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB);
    end
    
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
    
     %调试
     if AllConnected == 0
        fprintf('================AfterConnectALLAdj---PrintAM====================\n');
        am = AM_{ClusterIdx(1),ClusterIdx(2)};
        nodes = size(am,2);
        for i = 1:nodes
            for j = 1:nodes
                fprintf('am(%d,%d) = %d\t',i,j,am(i,j));
            end
            fprintf('\n');
        end
        fprintf('=========================\n');
    
    
        fprintf('================ReCheckConnectivity====================\n');
        RealConnected = CheckConnected(am,affectedNodes(1),nodeIdx);
        fprintf('RealConnected = %d\n',RealConnected);
        fprintf('=========================\n');
     end
    
    %节点离开
    [ClusterMatrix_,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_]...
    = NodeOut(nodeIdx,nodesOutNum,ClusterMatrix_,ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_);
    
    %%节点加入
    [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay] = ...
    NewNodeJoin(MaxLinkDistance,ClusterMatrix_,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexDelay_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,...
        EdgeBWDUB,VertexDelayDUB);
end

%节点离开
function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = NodeOut(nodeOutIdx,nodesOutNum,ClusterMatrix_,ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_)
 
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    
    %更新AM   
    am(:,nodeOutIdx) = [];
    am(nodeOutIdx,:) = [];
    AM_{ClusterIdx(1),ClusterIdx(2)} = am; 
    
    %更新ClusterMatrix
    Cluster(:,nodeOutIdx) = [];
    Cluster(1,:) = [1:nodesOutNum-1];
    ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)} = Cluster; 
    
    %更新Edge*
    edgeCost(:,nodeOutIdx) = [];
    edgeCost(nodeOutIdx,:) = [];
    EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;
    
    edgeDelay(:,nodeOutIdx) = [];
    edgeDelay(nodeOutIdx,:) = [];
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
    
    edgeBW(:,nodeOutIdx) = [];
    edgeBW(nodeOutIdx,:) = [];
    EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW;
    
    %更新Vertex*    
    vertexDelay(nodeOutIdx) = [];
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    
    ClusterMatrix = ClusterMatrix_;
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
end

%同一个簇节点重新加入
function [Connected,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = OldNodeJoin(MaxLinkDistance,affectedNodes,ClusterMatrix_,ClusterIdx,affectedNodeIdx,affectedNodePos,nodeIdx,AM_,...
        EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%% 从退出节点的相邻节点中选择距离最近的节点相连

    %选择与nodeIdx最近的邻结点相连
    sortedNodes = SortNodesByDistance(0,affectedNodePos,ClusterMatrix_,ClusterIdx,affectedNodes);
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    affectedNum = size(sortedNodes,2);
    Connected = 0;
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};  
    for i = 1:affectedNum
        linkIdx = sortedNodes(i);
        linkNode = Cluster(:,linkIdx); 
        if linkIdx ~= affectedNodeIdx && am(affectedNodeIdx,linkIdx) == 0
            Distance = ( (affectedNodePos(2)-linkNode(2))^2 +...
                        (affectedNodePos(3)-linkNode(3))^2)^0.5;
            if Distance > MaxLinkDistance
                break;
            end
            
            am(affectedNodeIdx,linkIdx) = 1;
            am(linkIdx,affectedNodeIdx) = 1;
            
            %调试
            fprintf('linkIdx = %d\n',linkIdx);
            nodes = size(am,2);
            for k = 1:nodes
                fprintf('am(%d,%d) = %d\n',affectedNodeIdx,k,am(affectedNodeIdx,k));
            end
            
            Connected = CheckConnected(am,affectedNodeIdx,nodeIdx);
            fprintf('Connected=%d\n',Connected);
            if Connected == 1
                break;
            else
                am(affectedNodeIdx,linkIdx) = 0;
                am(linkIdx,affectedNodeIdx) = 0;
            end
        end
    end
    
    %找到连接的结点使拓扑连通
    if Connected == 1       
        %更新AM
        AM_{ClusterIdx(1),ClusterIdx(2)} = am;

        %更新EdgeCost
        edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
        %TODO应有公式
        NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
        edgeCost(affectedNodeIdx,linkIdx) = NewEdgeCost;
        edgeCost(linkIdx,affectedNodeIdx) = edgeCost(affectedNodeIdx,linkIdx);
        EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;

        %更新EdgeDelay
        edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
        NewEdgeDelay = 0.5*Distance/100000;
        edgeDelay(affectedNodeIdx,linkIdx) = NewEdgeDelay;
        edgeDelay(linkIdx,affectedNodeIdx) = edgeDelay(affectedNodeIdx,linkIdx);
        EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;

        %更新EdgeBw
        edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
        NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
        edgeBW(affectedNodeIdx,linkIdx) = NewEdgeBW;
        edgeBW(linkIdx,affectedNodeIdx) = edgeBW(affectedNodeIdx,linkIdx);
        EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW;
        %更新VertexDelay
        vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
        NewVertexDelay = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
        vertexDelay(affectedNodeIdx) = NewVertexDelay;
        VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    end
    
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;

end

%另一个簇节点加入
function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = NewNodeJoin(MaxLinkDistance,ClusterMatrix_,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
       VertexDelay_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
    
    row = min(max(ceil(10 * nodePos(3)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * nodePos(2)/(BorderLength * ColBase)),1),ColCnt);

    %更新ClusterMatrix
    ClusterIn = ClusterMatrix_{row,col};
    nodesIn = size(ClusterIn,2);
    ClusterIn = [ClusterIn,[nodesIn+1,nodePos(2),nodePos(3),nodePos(4),nodePos(5)]'];
    ClusterMatrix_{row,col} = ClusterIn;
    ClusterMatrix = ClusterMatrix_;

    %选择距离最近的节点搭建链路(检查链路有效性）
    sortedNodes = SortNodesByDistance(1,nodePos,ClusterMatrix_,[row,col],[1:nodesIn]);
    linkIdx = sortedNodes(1);
    linkNode = ClusterIn(:,linkIdx);
    Distance = ( (nodePos(2)-linkNode(2))^2 +...
                    (nodePos(3)-linkNode(3))^2)^0.5;
     
    %更新AM
    amIn = AM_{row,col};
    amAddCol = zeros(nodesIn,1);
    amAddRow = zeros(1,nodesIn+1);
    if Distance < MaxLinkDistance
        amAddCol(linkIdx) = 1;
        amAddRow(linkIdx) = 1;
    else
        amAddCol(linkIdx) = 0;
        amAddRow(linkIdx) = 0;
    end
    amIn = [amIn,amAddCol];
    amIn(nodesIn+1,:) = amAddRow;
    AM_{row,col} = amIn;

    %更新EdgeCost
    edgeCostIn = EdgeCost_{row,col};
    edgeCostAddCol = Inf(nodesIn,1);
    edgeCostAddRow = Inf(1,nodesIn+1); 
    if Distance < MaxLinkDistance
        %TODO 应有公式
        NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
        edgeCostAddCol(linkIdx) = NewEdgeCost;
        edgeCostAddRow(linkIdx) = NewEdgeCost;
    else
        edgeCostAddCol(linkIdx) = inf;
        edgeCostAddRow(linkIdx) = inf;
    end
    edgeCostIn = [edgeCostIn,edgeCostAddCol];
    edgeCostIn(nodesIn+1,:) = edgeCostAddRow;
    EdgeCost_{row,col} = edgeCostIn;

    %更新EdgeDelay
    edgeDelayIn = EdgeDelay_{row,col};
    edgeDelayAddCol = Inf(nodesIn,1);
    edgeDelayAddRow = Inf(1,nodesIn+1);
    if Distance < MaxLinkDistance
        NewEdgeDelay = 0.5*Distance/100000;
        edgeDelayAddCol(linkIdx) = NewEdgeDelay;
        edgeDelayAddRow(linkIdx) = NewEdgeDelay;
    else
        edgeDelayAddCol(linkIdx) = inf;
        edgeDelayAddRow(linkIdx) = inf; 
    end   
    edgeDelayIn = [edgeDelayIn,edgeDelayAddCol];
    edgeDelayIn(nodesIn+1,:) = edgeDelayAddRow;
    EdgeDelay_{row,col} = edgeDelayIn;

    %更新EdgeBW
    edgeBWIn = EdgeBW_{row,col};
    edgeBWAddCol = Inf(nodesIn,1);
    edgeBWAddRow = Inf(1,nodesIn+1); 
    if Distance < MaxLinkDistance
        %TODO 应有公式
        NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
        edgeBWAddCol(linkIdx) = NewEdgeBW;
        edgeBWAddRow(linkIdx) = NewEdgeBW;
    else
        edgeBWAddCol(linkIdx) = inf;
        edgeBWAddRow(linkIdx) = inf;
    end
    edgeBWIn = [edgeBWIn,edgeBWAddCol];
    edgeBWIn(nodesIn+1,:) = edgeBWAddRow;
    EdgeBW_{row,col} = edgeBWIn;

    %更新Vertex*
    vertexDelayIn = VertexDelay_{row,col};
    if Distance < MaxLinkDistance
        %TODO 应有公式
        vertexDelayIn(nodesIn+1) = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
    else
        vertexDelayIn(nodesIn+1) = inf;
    end
    VertexDelay_{row,col} = vertexDelayIn;  
    
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
end

function [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay] = ...
        ConnectNodeIdxAllAdj(MaxLinkDistance,ClusterMatrix_,ClusterIdx,affectedNodes,AM_,EdgeCost_,EdgeDelay_,...
        EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%%将affectedNodes按其横坐标排序，然后首尾互连
    affectedNum = size(affectedNodes,2);
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    
    affectedNodesWithX = zeros(2,affectedNum);
    for i = 1:affectedNum
        affectedNodesWithX(:,i) = Cluster(1:2,affectedNodes(i));
    end
    X = affectedNodesWithX(2,:);
    sorted_X = sort(X);
    sortedAffectedNodes = zeros(1,affectedNum);
    for i = 1:affectedNum
        pos = find(X==sorted_X(i));
        sortedAffectedNodes(i) = affectedNodesWithX(1,pos);
    end
    
    %调试
    fprintf('==============ConnectNodeIdxAllAdj---sortedAffectedNodes===============\n');
    for i = 1:affectedNum
        fprintf('%d\t',sortedAffectedNodes(i));
    end
    fprintf('\n===================================================================\n');
    
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    node = sortedAffectedNodes(1);
    for i = 2:affectedNum
        %只处理不互连的结点(链路长度须有效）
        nextNode = sortedAffectedNodes(i);
        nodePos = Cluster(:,node);
        nextNodePos = Cluster(:,nextNode);
        Distance = ( (nodePos(2)-nextNodePos(2))^2 +...
                        (nodePos(3)-nextNodePos(3))^2)^0.5;
        if am(node,nextNode) == 0 && Distance < MaxLinkDistance
            %更新AM_
            am(node,nextNode) = 1;
            am(nextNode,node) = 1;

            %更新EdgeCost
            %TODO应有公式
            NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
            edgeCost(node,nextNode) = NewEdgeCost;
            edgeCost(nextNode,node) = edgeCost(node,nextNode);           

            %更新EdgeDelay
            NewEdgeDelay = 0.5*Distance/100000;
            edgeDelay(node,nextNode) = NewEdgeDelay;
            edgeDelay(nextNode,node) = edgeDelay(node,nextNode);

            %更新EdgeBw
            NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
            edgeBW(node,nextNode) = NewEdgeBW;
            edgeBW(nextNode,node) = edgeBW(node,nextNode);
                      
            %更新VertexDelay
            NewVertexDelay = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
            vertexDelay(node) = NewVertexDelay;        
        end
        node = nextNode;
    end
    
    %调试
    fprintf('================ConnectALLAdj---PrintAM====================\n');
    nodes = size(am,2);
    for i = 1:nodes
        for j = 1:nodes
            fprintf('am(%d,%d) = %d\t',i,j,am(i,j));
        end
        fprintf('\n');
    end
    fprintf('=========================\n');
    
    AM_{ClusterIdx(1),ClusterIdx(2)} = am;
    EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
    EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW; 
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;  
    
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
end


