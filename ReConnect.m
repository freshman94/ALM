function [IsClusterHead,AM,EdgeDelay,VertexDelay,LDT,VertexStability,VertexPriority]...
    = ReConnect(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,AM_,...
    EdgeDelay_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,VertexPriority_)
%%若节点脱离，尝试让节点重新连接
%%策略：从am(i,j) == 0的结点中寻找最近的结点，并与其连接
%%(1)若当前结点nodeIdx度超出，break
%%(2)若最近的结点未超出MaxLinkDistance,但度超出，则尝试与下一个最近的结点相连
%%(3)若最近结点已超出MaxLinkDistance，该结点无法重连，break
%%(4)若与最近的结点相连，拓扑通，重连成功，break
%%(5)若与最近结点相连，拓扑仍不通，则使最近的结点重复以上过程
%%(6)若所有结点均已重复以上过程，重连过程结束
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    isConnected = CheckConnected(am,nodeIdx,inf);
    if isConnected == 0
        fprintf('=============================ReConnect=========================\n');
        Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
        edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
        ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
        vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
        vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
        vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
        nodesNum = size(am,1);
        vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
        for i = 1:nodesNum
            nodeDegree = sum(am(nodeIdx,:));
            nodeMaxDegree = vertexMaxDegree(nodeIdx);
            if nodeDegree >= nodeMaxDegree
                break;
            end
            NotConnectedNodes = [];
            for k = 1:nodesNum
                if am(nodeIdx,k) == 0
                    NotConnectedNodes = [NotConnectedNodes,k];
                end
            end            
            
            nodePos = Cluster(:,nodeIdx);
            sortedNodes = SortNodesByDistance(0,nodePos,ClusterMatrix_,ClusterIdx,...
                NotConnectedNodes);
            Distance = inf;
                                  
            NotConnectedNodesNum = size(NotConnectedNodes,2);
            for p = 1:NotConnectedNodesNum
                linkIdx = sortedNodes(1,p);
                Distance = sortedNodes(2,p);
                linkNodeDegree = sum(am(linkIdx,:));
                linkNodeMaxDegree = vertexMaxDegree(linkIdx);
                if linkNodeDegree >= linkNodeMaxDegree
                    continue;
                else
                    break;
                end
            end
            
            if Distance > MaxLinkDistance
                break;
            else
                %添加新链路
                am(nodeIdx,linkIdx) = 1;
                am(linkIdx,nodeIdx) = 1;
                
                NewEdgeDelay = 0.5*Distance/100000;
                edgeDelay(nodeIdx,linkIdx) = NewEdgeDelay;
                edgeDelay(linkIdx,nodeIdx) = edgeDelay(nodeIdx,linkIdx);

                for k = 1:nodesNum
                    vertexDelay(k) = GetVertexDelay(k,am,edgeDelay);
                end
         
                ldt(nodeIdx,linkIdx) = GetLDT(MaxLinkDistance,Cluster(:,[nodeIdx,linkIdx]));
                ldt(linkIdx,nodeIdx) = ldt(nodeIdx,linkIdx);
                
                vertexStability(nodeIdx) = sum(ldt(nodeIdx,:));
                vertexStability(linkIdx) = sum(ldt(linkIdx,:));
                
                for k = 1:nodesNum
                    vertexPriority(i) = GetPriority(i,vertexStability,...
                        am,vertexDelay);
                end
                IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
                
                AM_{ClusterIdx(1),ClusterIdx(2)} = am;
                EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
                LDT_{ClusterIdx(1),ClusterIdx(2)} = ldt;
                VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
                VertexStability_{ClusterIdx(1),ClusterIdx(2)} = vertexStability;
                VertexPriority_{ClusterIdx(1),ClusterIdx(2)} = vertexPriority;
                
            end

            isConnected = CheckConnected(am,nodeIdx,inf);
            if isConnected == 1
                break;
            else
                nodeIdx = linkIdx;
            end
        end
    end
    
    IsClusterHead = IsClusterHead_;
    AM = AM_;
    EdgeDelay = EdgeDelay_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
end