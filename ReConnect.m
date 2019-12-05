function [IsClusterHead,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay,LDT,VertexStability,...
    VertexPriority]...
    = ReConnect(IsClusterHead_,MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,AM_,...
    EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,VertexMaxDegree_,LDT_,VertexStability_,...
    VertexPriority_,EdgeCostDUB,EdgeBWDUB)
%%���ڵ����룬�����ýڵ���������
%%���ԣ���am(i,j) == 0�Ľ����Ѱ������Ľ�㣬����������
%%(1)����ǰ���nodeIdx�ȳ�����break
%%(2)������Ľ��δ����MaxLinkDistance,���ȳ�������������һ������Ľ������
%%(3)���������ѳ���MaxLinkDistance���ý���޷�������break
%%(4)��������Ľ������������ͨ�������ɹ���break
%%(5)���������������������Բ�ͨ����ʹ����Ľ���ظ����Ϲ���
%%(6)�����н������ظ����Ϲ��̣��������̽���
    fprintf('=============================ReConnect=========================\n');

    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    isConnected = CheckConnected(am,nodeIdx,inf);
    if isConnected == 0
        Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
        edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
        edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
        ldt = LDT_{ClusterIdx(1),ClusterIdx(2)};
        vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
        vertexStability = VertexStability_{ClusterIdx(1),ClusterIdx(2)};
        vertexPriority = VertexPriority_{ClusterIdx(1),ClusterIdx(2)};
        edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
        nodesNum = size(am,1);
        vertexMaxDegree = VertexMaxDegree_{ClusterIdx(1),ClusterIdx(2)};
        for i = 1:nodesNum
            NotConnectedNodes = [];
            for k = 1:nodesNum
                if am(nodeIdx,k) == 0
                    NotConnectedNodes = [NotConnectedNodes,k];
                end
            end
            nodeDegree = sum(am(nodeIdx,:));
            nodeMaxDegree = vertexMaxDegree(nodeIdx);
            nodePos = Cluster(:,nodeIdx);
            sortedNodes = SortNodesByDistance(0,nodePos,ClusterMatrix_,ClusterIdx,NotConnectedNodes);
            linkIdx = sortedNodes(1);
            linkNodeDegree = sum(am(linkIdx,:));
            linkNodeMaxDegree = vertexMaxDegree(linkIdx);
            linkNodePos = Cluster(:,linkIdx);
            
            if nodeDegree >= nodeMaxDegree
                break;
            end
            
            NotConnectedNodesNum = size(NotConnectedNodes,2);
            for p = 2:NotConnectedNodesNum
                if linkNodeDegree >= linkNodeMaxDegree
                    linkIdx = sortedNodes(p);
                    linkNodeDegree = sum(am(linkIdx,:));
                    linkNodeMaxDegree = vertexMaxDegree(linkIdx);
                    linkNodePos = Cluster(:,linkIdx);
                    continue;
                else
                    break;
                end
            end
            
            Distance = ( (nodePos(2)-linkNodePos(2))^2 +...
                        (nodePos(3)-linkNodePos(3))^2)^0.5;
            if Distance > MaxLinkDistance
                break;
            else
                %�������·
                am(nodeIdx,linkIdx) = 1;
                am(linkIdx,nodeIdx) = 1;
                %TODOӦ�й�ʽ
                NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
                edgeCost(nodeIdx,linkIdx) = NewEdgeCost;
                edgeCost(linkIdx,nodeIdx) = edgeCost(nodeIdx,linkIdx);
                
                NewEdgeDelay = 0.5*Distance/100000;
                edgeDelay(nodeIdx,linkIdx) = NewEdgeDelay;
                edgeDelay(linkIdx,nodeIdx) = edgeDelay(nodeIdx,linkIdx);
                
                NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
                edgeBW(nodeIdx,linkIdx) = NewEdgeBW;
                edgeBW(linkIdx,nodeIdx) = edgeBW(nodeIdx,linkIdx);

                for k = 1:nodesNum
                    vertexDelay(k) = GetVertexDelay(k,am,edgeDelay);
                end
         
                ldt(nodeIdx,linkIdx) = GetLDT(MaxLinkDistance,Cluster(:,[nodeIdx,linkIdx]));
                ldt(linkIdx,nodeIdx) = ldt(nodeIdx,linkIdx);
                
                vertexStability(nodeIdx) = sum(ldt(nodeIdx,:));
                vertexStability(linkIdx) = sum(ldt(linkIdx,:));
                
                for k = 1:nodesNum
                    vertexPriority(i) = GetPriority(i,vertexStability,...
                        am,vertexMaxDegree,vertexDelay);
                end
                IsClusterHead_ = SetClusterHead(IsClusterHead_,ClusterIdx,vertexPriority);
                
                AM_{ClusterIdx(1),ClusterIdx(2)} = am;
                EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;
                EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
                EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW;  
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
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
    LDT = LDT_;
    VertexStability = VertexStability_;
    VertexPriority = VertexPriority_;
end