function [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = ReConnect(MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,AM_,EdgeCost_,EdgeDelay_,...
        EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%%若节点脱离，尝试让节点重新连接
%%策略：从am(i,j) == 0的结点中寻找最近的结点，并与其连接
%%(1)若最近结点已超出MaxLinkDistance，该结点无法重连，break
%%(2)若与最近的结点相连，拓扑通，重连成功，break
%%(3)若与最近结点相连，拓扑仍不通，则使最近的结点重复以上过程
%%(4)若所有结点均已重复以上过程，重连过程结束
    
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    isConnected = CheckConnected(am,nodeIdx,inf);
    if isConnected == 0
        Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
        nodesNum = size(am,1);
        for i = 1:nodesNum
            NotConnectedNodes = [];
            for k = 1:nodesNum
                if am(nodeIdx,k) == 0
                    NotConnectedNodes = [NotConnectedNodes,k];
                end
            end
            nodePos = Cluster(:,nodeIdx);
            sortedNodes = SortNodesByDistance(0,nodePos,ClusterMatrix_,ClusterIdx,NotConnectedNodes);
            linkIdx = sortedNodes(1);
            linkNode = Cluster(:,linkIdx);
            Distance = ( (nodePos(2)-linkNode(2))^2 +...
                        (nodePos(3)-linkNode(3))^2)^0.5;
            if Distance > MaxLinkDistance
                break;
            else
                %添加新链路
                am(nodeIdx,linkIdx) = 1;
                am(linkIdx,nodeIdx) = 1;
                AM_{ClusterIdx(1),ClusterIdx(2)} = am;
                
                edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
                %TODO应有公式
                NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
                edgeCost(nodeIdx,linkIdx) = NewEdgeCost;
                edgeCost(linkIdx,nodeIdx) = edgeCost(nodeIdx,linkIdx);
                
                edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
                NewEdgeDelay = 0.5*Distance/100000;
                edgeDelay(nodeIdx,linkIdx) = NewEdgeDelay;
                edgeDelay(linkIdx,nodeIdx) = edgeDelay(nodeIdx,linkIdx);
                
                edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
                NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
                edgeBW(nodeIdx,linkIdx) = NewEdgeBW;
                edgeBW(linkIdx,nodeIdx) = edgeBW(nodeIdx,linkIdx);
                
                vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
                NewVertexDelay = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
                vertexDelay(nodeIdx) = NewVertexDelay;
                
                AM_{ClusterIdx(1),ClusterIdx(2)} = am;
                EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;
                EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
                EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW; 
                VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;  
            end

            isConnected = CheckConnected(am,nodeIdx,inf);
            if isConnected == 1
                break;
            else
                nodeIdx = linkIdx;
            end
        end
    end
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
end