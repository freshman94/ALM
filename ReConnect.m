function [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = ReConnect(MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,AM_,EdgeCost_,EdgeDelay_,...
        EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%%���ڵ����룬�����ýڵ���������
%%���ԣ���am(i,j) == 0�Ľ����Ѱ������Ľ�㣬����������
%%(1)���������ѳ���MaxLinkDistance���ý���޷�������break
%%(2)��������Ľ������������ͨ�������ɹ���break
%%(3)���������������������Բ�ͨ����ʹ����Ľ���ظ����Ϲ���
%%(4)�����н������ظ����Ϲ��̣��������̽���
    
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
                %�������·
                am(nodeIdx,linkIdx) = 1;
                am(linkIdx,nodeIdx) = 1;
                AM_{ClusterIdx(1),ClusterIdx(2)} = am;
                
                edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
                %TODOӦ�й�ʽ
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