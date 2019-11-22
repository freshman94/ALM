function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss]...
    = VertexOutAndIn(ClusterMatrix_,ClusterIdx,nodeIdx,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexCost_,VertexDelay_,VertexPacketLoss_,BorderLength,RowBase,ColBase,EdgeCostDUB,...
        EdgeBWDUB,VertexCostDUB,VertexDelayDUB,VertexPacketLossDUB)
%%输入参数
%ClusterIdx: 1*2行向量，指明某个簇
%nodeIdx: 数值，指明簇中某个节点

    ClusterOut = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    
    %%节点离开，更新旧簇信息
    nodesOut = size(ClusterOut,2);
    %更新ClusterMatrix
    ClusterOut(:,nodeIdx) = [];
    ClusterOut(1,:) = [1:nodesOut-1];
    ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)} = ClusterOut;
    %更新AM
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    am(:,nodeIdx) = [];
    am(nodeIdx,:) = [];
    AM_{ClusterIdx(1),ClusterIdx(2)} = am;
    %更新Edge*
    edgeCostOut = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
    edgeCostOut(:,nodeIdx) = [];
    edgeCostOut(nodeIdx,:) = [];
    EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCostOut;
    
    edgeDelayOut = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelayOut(:,nodeIdx) = [];
    edgeDelayOut(nodeIdx,:) = [];
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelayOut;
    
    edgeBWOut = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
    edgeBWOut(:,nodeIdx) = [];
    edgeBWOut(nodeIdx,:) = [];
    EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBWOut;
    %更新Vertex*
    vertexCostOut = VertexCost_{ClusterIdx(1),ClusterIdx(2)};
    fprintf('nodes = %d, nodeIdx = %d\n',size(vertexCostOut,2),nodeIdx);
    vertexCostOut(nodeIdx) = [];
    VertexCost_{ClusterIdx(1),ClusterIdx(2)} = vertexCostOut;
    
    vertexDelayOut = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelayOut(nodeIdx) = [];
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelayOut;
    
    vertexPacketLossOut = VertexPacketLoss_{ClusterIdx(1),ClusterIdx(2)};
    vertexPacketLossOut(nodeIdx) = [];
    VertexPacketLoss_{ClusterIdx(1),ClusterIdx(2)} = vertexPacketLossOut;
    
    %%节点加入
    row = max(ceil(10 * nodePos(3)/(BorderLength*RowBase)),1);
    col = max(ceil(10 * nodePos(2)/(BorderLength * ColBase)),1);
    
    %更新ClusterMatrix
    ClusterIn = ClusterMatrix_{row,col};
    nodesIn = size(ClusterIn,2);
    ClusterIn = [ClusterIn,[nodesIn+1,nodePos(2),nodePos(3),nodePos(4),nodePos(5)]'];
    ClusterMatrix_{row,col} = ClusterIn;
    ClusterMatrix = ClusterMatrix_;
    
    %随机选择某个节点搭建链路
    linkIdx = randperm(nodesIn,1);
    linkNode = ClusterIn(:,linkIdx);
    
    %更新AM
    amIn = AM_{row,col};
    amAddCol = zeros(nodesIn,1);
    amAddRow = zeros(1,nodesIn+1);
    amAddCol(linkIdx) = 1;
    amAddRow(linkIdx) = 1;
    amIn = [amIn,amAddCol];
    amIn(nodesIn+1,:) = amAddRow;
    AM_{row,col} = amIn;
    AM = AM_;
    
    %更新EdgeCost
    edgeCostIn = EdgeCost_{row,col};
    %TODO 应有公式
    NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
    edgeCostAddCol = Inf(nodesIn,1);
    edgeCostAddRow = Inf(1,nodesIn+1);
    edgeCostAddCol(linkIdx) = NewEdgeCost;
    edgeCostAddRow(linkIdx) = NewEdgeCost;
    edgeCostIn = [edgeCostIn,edgeCostAddCol];
    edgeCostIn(nodesIn+1,:) = edgeCostAddRow;
    EdgeCost_{row,col} = edgeCostIn;
    EdgeCost = EdgeCost_;
    
    %更新EdgeDelay
    edgeDelayIn = EdgeDelay_{row,col};
    Distance = ( (nodePos(2)-linkNode(2))^2 +...
                    (nodePos(3)-linkNode(3))^2)^0.5;
    NewEdgeDelay = 0.5*Distance/100000;
    edgeDelayAddCol = Inf(nodesIn,1);
    edgeDelayAddRow = Inf(1,nodesIn+1);
    edgeDelayAddCol(linkIdx) = NewEdgeDelay;
    edgeDelayAddRow(linkIdx) = NewEdgeDelay;
    edgeDelayIn = [edgeDelayIn,edgeDelayAddCol];
    edgeDelayIn(nodesIn+1,:) = edgeDelayAddRow;
    EdgeDelay_{row,col} = edgeDelayIn;
    EdgeDelay = EdgeDelay_;
    
    %更新EdgeBW
    edgeBWIn = EdgeBW_{row,col};
    %TODO 应有公式
    NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
    edgeBWAddCol = Inf(nodesIn,1);
    edgeBWAddRow = Inf(1,nodesIn+1);
    edgeBWAddCol(linkIdx) = NewEdgeBW;
    edgeBWAddRow(linkIdx) = NewEdgeBW;
    edgeBWIn = [edgeBWIn,edgeBWAddCol];
    edgeBWIn(nodesIn+1,:) = edgeBWAddRow;
    EdgeBW_{row,col} = edgeBWIn;
    EdgeBW = EdgeBW_;
    
    %更新Vertex*
    vertexCostIn = VertexCost_{row,col};
    %TODO 应有公式
    vertexCostIn(nodesIn+1) = VertexCostDUB(1)+(VertexCostDUB(2)-VertexCostDUB(1))*rand;
    VertexCost_{row,col} = vertexCostIn;
    VertexCost = VertexCost_;
    
    vertexDelayIn = VertexDelay_{row,col};
    %TODO 应有公式
    vertexDelayIn(nodesIn+1) = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
    VertexDelay_{row,col} = vertexDelayIn;
    VertexDelay = VertexDelay_;
    
    vertexPacktLossIn = VertexPacketLoss_{row,col};
    %TODO 应有公式
    vertexPacktLossIn(nodesIn+1) = VertexPacketLossDUB(1)+(VertexPacketLossDUB(2)-VertexPacketLossDUB(1))*rand;
    VertexPacketLoss_{row,col} = vertexPacktLossIn;
    VertexPacketLoss = VertexPacketLoss_;
    
    %调试代码
    [row,col] = size(EdgeCost);
    for i = 1:row
        for j = 1:col
            nodes = size(EdgeCost{i,j},2);
            fprintf('In nodes = %d\n',nodes);
        end
    end
    
end