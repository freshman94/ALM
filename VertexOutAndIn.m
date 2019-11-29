function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = VertexOutAndIn(MaxLinkDistance,ClusterMatrix_,ClusterIdx,nodeIdx,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexDelay_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%%�������
%ClusterIdx: 1*2��������ָ��ĳ����
%nodeIdx: ��ֵ��ָ������ĳ���ڵ�

    ClusterOut = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    nodesOutNum = size(ClusterOut,2);
    %ͳ����Ӱ��ڵ�
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    %�ж���Ӱ��ڵ��Ƿ���Ҫ���¼����
    affectedNodes = [];
    for i = 1:nodesOutNum
        if am(nodeIdx,i) == 1 
            affectedNodes = [affectedNodes,i];
        end
    end  
    
    affectedNum = size(affectedNodes,2);
    
    %����
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
        %���ڵ�Ķ�Ϊ1�������˳���������������Ӱ��
        isConnected = CheckConnected(am,affectedNodes(i),nodeIdx);
        %���Դ���
        fprintf('affectedNode = %d,isConnected=%d\n',affectedNodes(i),isConnected);
        %�˱�������ָʾ����Ľڵ��Ƿ����¼���ɹ�
        %�ýڵ����룬���¼���
        if isConnected == 0
            affectedNodePos = ClusterOut(:,affectedNodes(i));
            [AllConnected,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_] = ...
            OldNodeJoin(MaxLinkDistance,affectedNodes,ClusterMatrix_,ClusterIdx,affectedNodes(i),affectedNodePos,nodeIdx,AM_,...
            EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB);
            %������ͨ���ؽ��ɹ�
%             if AllConnected == 1
%                 break;
%             end
        end
    end
    
    fprintf('AllConnected = %d\n',AllConnected);
    
    %����Ľڵ�û�м���ɹ�����nodeIdx�������ڽڵ������������¹���������ͨ��
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
    
    %�ڵ��뿪
    [ClusterMatrix_,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_]...
    = NodeOut(nodeIdx,nodesOutNum,ClusterMatrix_,ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_);
    
    %%�ڵ����
    [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay] = ...
    NewNodeJoin(MaxLinkDistance,ClusterMatrix_,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
        VertexDelay_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,...
        EdgeBWDUB,VertexDelayDUB);
end

%�ڵ��뿪
function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = NodeOut(nodeOutIdx,nodesOutNum,ClusterMatrix_,ClusterIdx,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_)
 
    Cluster = ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)};
    am = AM_{ClusterIdx(1),ClusterIdx(2)};
    edgeCost = EdgeCost_{ClusterIdx(1),ClusterIdx(2)};
    edgeDelay = EdgeDelay_{ClusterIdx(1),ClusterIdx(2)};
    edgeBW = EdgeBW_{ClusterIdx(1),ClusterIdx(2)};
    vertexDelay = VertexDelay_{ClusterIdx(1),ClusterIdx(2)};
    
    %����AM   
    am(:,nodeOutIdx) = [];
    am(nodeOutIdx,:) = [];
    AM_{ClusterIdx(1),ClusterIdx(2)} = am; 
    
    %����ClusterMatrix
    Cluster(:,nodeOutIdx) = [];
    Cluster(1,:) = [1:nodesOutNum-1];
    ClusterMatrix_{ClusterIdx(1),ClusterIdx(2)} = Cluster; 
    
    %����Edge*
    edgeCost(:,nodeOutIdx) = [];
    edgeCost(nodeOutIdx,:) = [];
    EdgeCost_{ClusterIdx(1),ClusterIdx(2)} = edgeCost;
    
    edgeDelay(:,nodeOutIdx) = [];
    edgeDelay(nodeOutIdx,:) = [];
    EdgeDelay_{ClusterIdx(1),ClusterIdx(2)} = edgeDelay;
    
    edgeBW(:,nodeOutIdx) = [];
    edgeBW(nodeOutIdx,:) = [];
    EdgeBW_{ClusterIdx(1),ClusterIdx(2)} = edgeBW;
    
    %����Vertex*    
    vertexDelay(nodeOutIdx) = [];
    VertexDelay_{ClusterIdx(1),ClusterIdx(2)} = vertexDelay;
    
    ClusterMatrix = ClusterMatrix_;
    AM = AM_;
    EdgeCost = EdgeCost_;
    EdgeDelay = EdgeDelay_;
    EdgeBW = EdgeBW_;
    VertexDelay = VertexDelay_;
end

%ͬһ���ؽڵ����¼���
function [Connected,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = OldNodeJoin(MaxLinkDistance,affectedNodes,ClusterMatrix_,ClusterIdx,affectedNodeIdx,affectedNodePos,nodeIdx,AM_,...
        EdgeCost_,EdgeDelay_,EdgeBW_,VertexDelay_,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
%% ���˳��ڵ�����ڽڵ���ѡ���������Ľڵ�����

    %ѡ����nodeIdx������ڽ������
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

%��һ���ؽڵ����
function [ClusterMatrix,AM,EdgeCost,EdgeDelay,EdgeBW,VertexDelay]...
    = NewNodeJoin(MaxLinkDistance,ClusterMatrix_,nodePos,AM_,EdgeCost_,EdgeDelay_,EdgeBW_,...
       VertexDelay_,BorderLength,RowBase,ColBase,RowCnt,ColCnt,EdgeCostDUB,EdgeBWDUB,VertexDelayDUB)
    
    row = min(max(ceil(10 * nodePos(3)/(BorderLength*RowBase)),1),RowCnt);
    col = min(max(ceil(10 * nodePos(2)/(BorderLength * ColBase)),1),ColCnt);

    %����ClusterMatrix
    ClusterIn = ClusterMatrix_{row,col};
    nodesIn = size(ClusterIn,2);
    ClusterIn = [ClusterIn,[nodesIn+1,nodePos(2),nodePos(3),nodePos(4),nodePos(5)]'];
    ClusterMatrix_{row,col} = ClusterIn;
    ClusterMatrix = ClusterMatrix_;

    %ѡ���������Ľڵ���·(�����·��Ч�ԣ�
    sortedNodes = SortNodesByDistance(1,nodePos,ClusterMatrix_,[row,col],[1:nodesIn]);
    linkIdx = sortedNodes(1);
    linkNode = ClusterIn(:,linkIdx);
    Distance = ( (nodePos(2)-linkNode(2))^2 +...
                    (nodePos(3)-linkNode(3))^2)^0.5;
     
    %����AM
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

    %����EdgeCost
    edgeCostIn = EdgeCost_{row,col};
    edgeCostAddCol = Inf(nodesIn,1);
    edgeCostAddRow = Inf(1,nodesIn+1); 
    if Distance < MaxLinkDistance
        %TODO Ӧ�й�ʽ
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

    %����EdgeDelay
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

    %����EdgeBW
    edgeBWIn = EdgeBW_{row,col};
    edgeBWAddCol = Inf(nodesIn,1);
    edgeBWAddRow = Inf(1,nodesIn+1); 
    if Distance < MaxLinkDistance
        %TODO Ӧ�й�ʽ
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

    %����Vertex*
    vertexDelayIn = VertexDelay_{row,col};
    if Distance < MaxLinkDistance
        %TODO Ӧ�й�ʽ
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
        %ֻ���������Ľ��(��·��������Ч��
        nextNode = sortedAffectedNodes(i);
        nodePos = Cluster(:,node);
        nextNodePos = Cluster(:,nextNode);
        Distance = ( (nodePos(2)-nextNodePos(2))^2 +...
                        (nodePos(3)-nextNodePos(3))^2)^0.5;
        if am(node,nextNode) == 0 && Distance < MaxLinkDistance
            %����AM_
            am(node,nextNode) = 1;
            am(nextNode,node) = 1;

            %����EdgeCost
            %TODOӦ�й�ʽ
            NewEdgeCost = EdgeCostDUB(1)+(EdgeCostDUB(2)-EdgeCostDUB(1))*rand;
            edgeCost(node,nextNode) = NewEdgeCost;
            edgeCost(nextNode,node) = edgeCost(node,nextNode);           

            %����EdgeDelay
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


