function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss]...
    = VertexOutAndIn(ClusterMatrix_,ClusterIdx,nodeIdx,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexCost_,VertexDelay_,VertexPacketLoss_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,...
        EdgeBWDUB,VertexCostDUB,VertexDelayDUB,VertexPacketLossDUB)
%%�������
%ClusterIdx: 1*2��������ָ��ĳ����
%nodeIdx: ��ֵ��ָ������ĳ���ڵ�

    ClusterOut = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    
    %%�ڵ��뿪�����¾ɴ���Ϣ
    nodesOut = size(ClusterOut,2);
    %ͳ����Ӱ��ڵ�
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    %�ж���Ӱ��ڵ��Ƿ���Ҫ���¼����
    affectedNodes = [];
    for i = 1:nodesOut
        if am(nodeIdx,i) == 1 
            affectedNodes = [affectedNodes,i];
        end
    end  
    
    affectedNum = size(affectedNodes,2);
    
    %����
    fprintf('clusterIdx=[%d,%d],nodeIdx = %d,affectedNum = %d\n',ClusterIdx(1),ClusterIdx(2),nodeIdx,affectedNum);
    
    fprintf('================print am==========\n');
    for i = 1:nodesOut
        for j = 1:nodesOut
            fprintf('am(%d,%d) = %d\t',i,j,am(i,j));
        end
        fprintf('\n');
    end
    fprintf('=========================\n');
    
    am(:,nodeIdx) = zeros(nodesOut,1);
    am(nodeIdx,:) = zeros(1,nodesOut);
    AM_{ClusterIdx(1),ClusterIdx(2)} = am;
    
    AllConnected = 1;
    for i = 1:affectedNum
        am = AM_{ClusterIdx(1),ClusterIdx(2)};
        %���ڵ�Ķ�Ϊ1�������˳���������������Ӱ��
        isConnected = CheckConnected(am,affectedNodes(i),nodeIdx);
        %���Դ���
        fprintf('affectedNode = %d,isConnected=%d\n',affectedNodes(i),isConnected);
        %�˱�������ָʾ����Ľڵ��Ƿ����¼���ɹ�
        %�ýڵ����룬���¼���
        if isConnected == 0
            affectedNodePos = ClusterOut(:,affectedNodes(i));
            [AllConnected,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_] = ...
            OldNodeJoin(affectedNodes,ClusterMatrix_,ClusterIdx,affectedNodes(i),affectedNodePos,nodeIdx,AM_,...
            EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB);           
        end
    end
    
    fprintf('AllConnected = %d\n',AllConnected);
    
    %����Ľڵ�û�м���ɹ�����nodeIdx�������ڽڵ������������¹���������ͨ��
    if AllConnected == 0
        [AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_] = ...
        ConnectNodeIdxAllAdj(ClusterMatrix_,ClusterIdx,affectedNodes,AM_,...
            EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB);
    end
    
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
    
     %����
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
    %����AM   
    am(:,nodeIdx) = [];
    am(nodeIdx,:) = [];
    AM_{ClusterIdx(1),ClusterIdx(2)} = am; 
    
    %����ClusterMatrix
    ClusterOut(:,nodeIdx) = [];
    ClusterOut(1,:) = [1:nodesOut-1];
    ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)} = ClusterOut; 
    
    %����Edge*
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
    %����Vertex*
    vertexCostOut = VertexCost_{ClusterIdx(1),ClusterIdx(2)};
    vertexCostOut(nodeIdx) = [];
    VertexCost_{ClusterIdx(1),ClusterIdx(2)} = vertexCostOut;
    
    vertexDelayOut = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelayOut(nodeIdx) = [];
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelayOut;
    
    vertexPacketLossOut = VertexPacketLoss_{ClusterIdx(1),ClusterIdx(2)};
    vertexPacketLossOut(nodeIdx) = [];
    VertexPacketLoss_{ClusterIdx(1),ClusterIdx(2)} = vertexPacketLossOut;
    
    %%�ڵ����
    [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss] = ...
    NewNodeJoin(ClusterMatrix_,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexCost_,VertexDelay_,VertexPacketLoss_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,...
        EdgeBWDUB,VertexCostDUB,VertexDelayDUB,VertexPacketLossDUB);
end


%ͬһ���ؽڵ����¼���
function [Connected,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = OldNodeJoin(affectedNodes,ClusterMatrix_,ClusterIdx,affectedNodeIdx,affectedNodePos,nodeIdx,AM_,...
        EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%% ���˳��ڵ�����ڽڵ���ѡ���������Ľڵ�����

    %ѡ����nodeIdx������ڽ������
    sortedNodes = SortNodesByDistance(0,affectedNodePos,ClusterMatrix_,ClusterIdx,affectedNodes);
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    affectedNum = size(sortedNodes,2);
    Connected = 0;
    for i = 1:affectedNum
        linkIdx = sortedNodes(i);
        if linkIdx ~= affectedNodeIdx && am(affectedNodeIdx,linkIdx) == 0
            am(affectedNodeIdx,linkIdx) = 1;
            am(linkIdx,affectedNodeIdx) = 1;
            
            %����
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
    
    %�ҵ����ӵĽ��ʹ������ͨ
    if Connected == 1       
        Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};   
        linkNode = Cluster(:,linkIdx);    

        %����AM
        AM_{ClusterIdx(1),ClusterIdx(2)} = am;

        %����EdgeCost
        edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
        %TODOӦ�й�ʽ
        NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
        edgeCost(affectedNodeIdx,linkIdx) = NewEdgeCost;
        edgeCost(linkIdx,affectedNodeIdx) = edgeCost(affectedNodeIdx,linkIdx);
        EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;

        %����EdgeDelay
        edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
        Distance = ( (affectedNodePos(2)-linkNode(2))^2 +...
                        (affectedNodePos(3)-linkNode(3))^2)^0.5;
        NewEdgeDelay = 0.5*Distance/100000;
        edgeDelay(affectedNodeIdx,linkIdx) = NewEdgeDelay;
        edgeDelay(linkIdx,affectedNodeIdx) = edgeDelay(affectedNodeIdx,linkIdx);
        EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;

        %����EdgeBw
        edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
        NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
        edgeBW(affectedNodeIdx,linkIdx) = NewEdgeBW;
        edgeBW(linkIdx,affectedNodeIdx) = edgeBW(affectedNodeIdx,linkIdx);
        EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW;
        %����VertexDelay
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

function [AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay] = ...
        ConnectNodeIdxAllAdj(ClusterMatrix_,ClusterIdx,affectedNodes,AM_,EdgeCost_,EdgeDelay_,...
        EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%%��affectedNodes�������������Ȼ����β����
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
    
    %����
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
        %ֻ���������Ľ��
        nextNode = sortedAffectedNodes(i);
        if am(node,nextNode) == 0
            %����AM_
            am(node,nextNode) = 1;
            am(nextNode,node) = 1;

            %����EdgeCost
            %TODOӦ�й�ʽ
            NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
            edgeCost(node,nextNode) = NewEdgeCost;
            edgeCost(nextNode,node) = edgeCost(node,nextNode);           

            %����EdgeDelay
            nodePos = Cluster(:,node);
            nextNodePos = Cluster(:,nextNode);
            Distance = ( (nodePos(2)-nextNodePos(2))^2 +...
                            (nodePos(3)-nextNodePos(3))^2)^0.5;
            NewEdgeDelay = 0.5*Distance/100000;
            edgeDelay(node,nextNode) = NewEdgeDelay;
            edgeDelay(nextNode,node) = edgeDelay(node,nextNode);

            %����EdgeBw
            NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
            edgeBW(node,nextNode) = NewEdgeBW;
            edgeBW(nextNode,node) = edgeBW(node,nextNode);
                      
            %����VertexDelay
            NewVertexDelay = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
            vertexDelay(node) = NewVertexDelay;        
        end
        node = nextNode;
    end
    
    %����
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

%�����ͨ��
function [isConnected] = CheckConnected(am,idx,nodeIdx)
    nodes = size(am,1);
    visit = zeros(1,nodes);
    fprintf('================================================\n');
    fprintf('DFS visit: ');
    visit = DFS(am,idx,visit,nodes);
    fprintf('\n');
    %����
    fprintf('visit: ');
    for i = 1:nodes
        fprintf('%d = %d\t',i,visit(i));
    end
    fprintf('\n================================================\n');
    
    for i = 1:nodes
        if  i ~= nodeIdx && visit(i) == 0
            isConnected = 0;
            return;
        end
    end
    isConnected = 1;
end

function [visit] = DFS (am,idx,visit_,nodes)
    visit_(idx) = 1;
    
    for i = 1:nodes
        %�ҵ����ӵĽڵ㣬�����ýڵ�
        if am(idx,i) ~= 0 && visit_(i) == 0 
            visit_(i) = 1;
            fprintf('%d = %d \t',i,visit_(i));
            visit_ = DFS(am,i,visit_,nodes);
        end      
    end
    visit = visit_;
end

%��һ���ؽڵ����
function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexCost,VertexDelay,VertexPacketLoss]...
    = NewNodeJoin(ClusterMatrix_,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexCost_,VertexDelay_,VertexPacketLoss_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,...
        EdgeBWDUB,VertexCostDUB,VertexDelayDUB,VertexPacketLossDUB)
    
    row = min(max(ceil(10 * nodePos(3)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * nodePos(2)/(BorderLength * ColBase)),1),ColCnt);

    %����ClusterMatrix
    ClusterIn = ClusterMatrix_{row,col};
    nodesIn = size(ClusterIn,2);
    ClusterIn = [ClusterIn,[nodesIn+1,nodePos(2),nodePos(3),nodePos(4),nodePos(5)]'];
    ClusterMatrix_{row,col} = ClusterIn;
    ClusterMatrix = ClusterMatrix_;

    %ѡ���������Ľڵ���·
    sortedNodes = SortNodesByDistance(1,nodePos,ClusterMatrix_,[row,col],[1:nodesIn]);
    linkIdx = sortedNodes(1);
    linkNode = ClusterIn(:,linkIdx);

    %����AM
    amIn = AM_{row,col};
    amAddCol = zeros(nodesIn,1);
    amAddRow = zeros(1,nodesIn+1);
    amAddCol(linkIdx) = 1;
    amAddRow(linkIdx) = 1;
    amIn = [amIn,amAddCol];
    amIn(nodesIn+1,:) = amAddRow;
    AM_{row,col} = amIn;
    AM = AM_;

    %����EdgeCost
    edgeCostIn = EdgeCost_{row,col};
    %TODO Ӧ�й�ʽ
    NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
    edgeCostAddCol = Inf(nodesIn,1);
    edgeCostAddRow = Inf(1,nodesIn+1);
    edgeCostAddCol(linkIdx) = NewEdgeCost;
    edgeCostAddRow(linkIdx) = NewEdgeCost;
    edgeCostIn = [edgeCostIn,edgeCostAddCol];
    edgeCostIn(nodesIn+1,:) = edgeCostAddRow;
    EdgeCost_{row,col} = edgeCostIn;
    EdgeCost = EdgeCost_;

    %����EdgeDelay
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

    %����EdgeBW
    edgeBWIn = EdgeBW_{row,col};
    %TODO Ӧ�й�ʽ
    NewEdgeBW = EdgeBWDUB(1)+(EdgeBWDUB(2)-EdgeBWDUB(1))*rand;
    edgeBWAddCol = Inf(nodesIn,1);
    edgeBWAddRow = Inf(1,nodesIn+1);
    edgeBWAddCol(linkIdx) = NewEdgeBW;
    edgeBWAddRow(linkIdx) = NewEdgeBW;
    edgeBWIn = [edgeBWIn,edgeBWAddCol];
    edgeBWIn(nodesIn+1,:) = edgeBWAddRow;
    EdgeBW_{row,col} = edgeBWIn;
    EdgeBW = EdgeBW_;

    %����Vertex*
    vertexCostIn = VertexCost_{row,col};
    %TODO Ӧ�й�ʽ
    vertexCostIn(nodesIn+1) = VertexCostDUB(1)+(VertexCostDUB(2)-VertexCostDUB(1))*rand;
    VertexCost_{row,col} = vertexCostIn;
    VertexCost = VertexCost_;

    vertexDelayIn = VertexDelay_{row,col};
    %TODO Ӧ�й�ʽ
    vertexDelayIn(nodesIn+1) = VertexDelayDUB(1)+(VertexDelayDUB(2)-VertexDelayDUB(1))*rand;
    VertexDelay_{row,col} = vertexDelayIn;
    VertexDelay = VertexDelay_;

    vertexPacktLossIn = VertexPacketLoss_{row,col};
    %TODO Ӧ�й�ʽ
    vertexPacktLossIn(nodesIn+1) = VertexPacketLossDUB(1)+(VertexPacketLossDUB(2)-VertexPacketLossDUB(1))*rand;
    VertexPacketLoss_{row,col} = vertexPacktLossIn;
    VertexPacketLoss = VertexPacketLoss_;   
end

function [sortedNodes] = SortNodesByDistance(NotIn,targetNodePos,ClusterMatrix_,ClusterIdx,nodes)
%% ��nodes�еĽڵ㰴����targeNode�ľ�����������
%NotIn����ָʾtargeNodePos�Ƿ��п��ܱ�nodes������1��ʾû�п���
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    nodeNum = size(nodes,2);
    nodesWithD = [];
    for i = 1:nodeNum
        if NotIn == 1 || targetNodePos(1) ~= nodes(i)
            distance = ( (targetNodePos(2)-Cluster(2,nodes(i)))^2 +...
                        (targetNodePos(3)-Cluster(3,nodes(i)))^2)^0.5;
            nodesWithD = [nodesWithD,[nodes(i),distance]'];
        end
    end
    sortedNodes = [];
    %nodes���ܰ���targeNode�������Ҫ����size
    nodeNum = size(nodesWithD,2);
    if nodeNum > 0
        distanceRow = nodesWithD(2,:);
        sortedNodesWithD = sort(distanceRow);
        sortedNodes = zeros(1,nodeNum);
        for i = 1:nodeNum
            pos = find(distanceRow==sortedNodesWithD(i));
            sortedNodes(i) = nodesWithD(1,pos);
        end
    end
        
    %����
    fprintf('============================SortNodesByDistance======================\n');
    for i = 1:nodeNum
        fprintf('%d\t',sortedNodes(i));
    end
    fprintf('\n===================================================================\n');
end
